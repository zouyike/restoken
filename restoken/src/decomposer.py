"""Inverse of the assembler: a peptide molecule -> ResToken sequence.

The assembler (cyclic_assembler.py) goes block IDs -> SMILES. This module
goes the other way: given a peptide SMILES (e.g. an approved macrocyclic
drug), it

  1. locates the peptide-backbone amide bonds,
  2. fragments the molecule into ordered residues,
  3. re-caps each residue as the Ac-/NHMe-capped form used by the library
     (so a fragment is directly comparable to a block's aa_smiles),
  4. matches each residue to its nearest library block by chirality-aware
     Morgan-fingerprint Tanimoto.

Backbone detection is amide-centric and chain-length aware: a backbone amide
C(=O)-N is one whose N reaches another backbone carbonyl through a short
all-carbon chain (1 carbon = alpha residue, 2 = beta, 3 = gamma). This both
handles non-alpha residues and rejects side-chain amides (Asn/Gln, linkers),
whose carbonyls are farther from the backbone nitrogen.

Side-chain crosslinks (e.g. the aromatic bridge that makes MK-0616 bicyclic)
leave two residues fused in one fragment; they are detected, broken, and
reported — head-to-tail ResToken cannot encode them, so they are an explicit,
quantified limitation rather than a silent error.

Validate the method with `roundtrip`, which assembles a known ResToken
sequence and checks the decomposer recovers the original block IDs.
"""

from collections import deque
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

from restoken.src.library import BlockLibrary


class PeptideDecomposer:
    def __init__(self, library: Optional[BlockLibrary] = None, radius: int = 2,
                 nbits: int = 2048, max_backbone_len: int = 3):
        self.lib = library or BlockLibrary()
        self.radius = radius
        self.nbits = nbits
        self.max_backbone_len = max_backbone_len
        # Sorted iteration makes tie-breaking among equal-Tanimoto blocks
        # (near-duplicate isomers) deterministic across runs.
        self._block_fp = {}
        for bid in sorted(self.lib.all_ids):
            m = Chem.MolFromSmiles(self.lib[bid].aa_smiles)
            if m is not None:
                self._block_fp[bid] = AllChem.GetMorganFingerprintAsBitVect(
                    m, radius, nBits=nbits, useChirality=True)

    # ── backbone analysis ────────────────────────────────────────────

    @staticmethod
    def _is_carbonyl_c(mol, idx):
        a = mol.GetAtomWithIdx(idx)
        if a.GetAtomicNum() != 6:
            return False
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                if mol.GetBondBetweenAtoms(idx, nb.GetIdx()).GetBondTypeAsDouble() == 2.0:
                    return True
        return False

    @staticmethod
    def _shares_small_ring(mol, a, b, max_size=9):
        """True if a and b share a ring no larger than max_size.

        Side-chain lactam rings are small (5-7); the head-to-tail backbone
        macrocycle is large (>=12 even for a 4-residue alpha peptide), so a
        small shared ring marks an intra-residue amide, not a backbone bond.
        """
        for ring in mol.GetRingInfo().AtomRings():
            if len(ring) <= max_size and a in ring and b in ring:
                return True
        return False

    def _amide_bonds(self, mol):
        """All C(=O)-N bonds as (carbonyl_c, n), excluding intra-ring lactams.

        An amide whose carbonyl and nitrogen share a ring is a side-chain
        lactam (e.g. the polyketo-proline blocks in this library), never a
        backbone peptide bond — the backbone carbonyl of a proline-like
        residue is exocyclic to its ring.
        """
        out = []
        seen = set()
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            for co, n in ((i, j), (j, i)):
                if (mol.GetAtomWithIdx(n).GetAtomicNum() == 7
                        and self._is_carbonyl_c(mol, co)
                        and (co, n) not in seen
                        and not self._shares_small_ring(mol, co, n)):
                    out.append((co, n))
                    seen.add((co, n))
        return out

    def _nearest_backbone_carbonyl(self, mol, n_idx, exclude_co, carbonyls):
        """Shortest all-carbon path from N to a carbonyl carbon (!= exclude_co)
        within max_backbone_len carbons. Returns that carbonyl idx or None."""
        # BFS over carbon atoms starting from neighbors of N
        start = n_idx
        dist = {start: 0}
        dq = deque([start])
        while dq:
            cur = dq.popleft()
            if dist[cur] >= self.max_backbone_len + 1:
                continue
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                j = nb.GetIdx()
                if j in dist:
                    continue
                if cur == start:
                    # first hop: step onto a carbon (the alpha/beta backbone C)
                    if nb.GetAtomicNum() != 6:
                        continue
                    dist[j] = 1
                    if j != exclude_co and j in carbonyls:
                        return j
                    dq.append(j)
                else:
                    # continue only through carbons; a carbonyl carbon is the target
                    if j in carbonyls and j != exclude_co:
                        return j
                    if nb.GetAtomicNum() == 6 and not self._is_carbonyl_c(mol, j):
                        dist[j] = dist[cur] + 1
                        dq.append(j)
        return None

    def _backbone_amides(self, mol):
        """Backbone amide bonds (co, n): N reaches another backbone carbonyl
        through a short carbon chain. Returns list and the residue adjacency
        map (carbonyl_of_residue -> carbonyl_of_next_residue)."""
        amides = self._amide_bonds(mol)
        carbonyls = {co for (co, _n) in amides}
        backbone = []
        for (co, n) in amides:
            co2 = self._nearest_backbone_carbonyl(mol, n, co, carbonyls)
            if co2 is not None:
                backbone.append((co, n, co2))
        return backbone

    # ── fragmentation / capping ──────────────────────────────────────

    def _capped_fragments(self, mol, backbone, crosslinks):
        """Break backbone + crosslink bonds on the intact molecule, adding
        Ac-/NHMe caps, then fragment. Operating on the perceived molecule
        preserves aromaticity. Returns (frag_idx_tuples, frag_smiles_list)."""
        rw = Chem.RWMol(mol)

        def add_nhme(c_idx):
            n = rw.AddAtom(Chem.Atom(7)); me = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(c_idx, n, Chem.BondType.SINGLE)
            rw.AddBond(n, me, Chem.BondType.SINGLE)

        def add_acetyl(n_idx):
            c = rw.AddAtom(Chem.Atom(6)); o = rw.AddAtom(Chem.Atom(8))
            me = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(n_idx, c, Chem.BondType.SINGLE)
            rw.AddBond(c, o, Chem.BondType.DOUBLE)
            rw.AddBond(c, me, Chem.BondType.SINGLE)

        for (co, n, _co2) in backbone:
            if rw.GetBondBetweenAtoms(co, n) is not None:
                rw.RemoveBond(co, n)
            add_nhme(co)
            add_acetyl(n)
        for (u, v) in crosslinks:
            if rw.GetBondBetweenAtoms(u, v) is not None:
                rw.RemoveBond(u, v)

        m = rw.GetMol()
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return [], []
        frag_idx = Chem.GetMolFrags(m, asMols=False, sanitizeFrags=False)
        frag_mols = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
        return frag_idx, [Chem.MolToSmiles(fm) for fm in frag_mols]

    # ── crosslink handling ───────────────────────────────────────────

    def _detect_crosslinks(self, mol, backbone):
        """After breaking backbone amides, any fragment holding >1 backbone
        carbonyl was fused by a side-chain crosslink. Break a non-ring bond on
        the path between the two carbonyls. Returns list of broken (u,v)."""
        bb_break = {tuple(sorted((co, n))) for (co, n, _c) in backbone}
        crosslinks = []
        for _ in range(20):
            comps = _components(mol, bb_break | {tuple(sorted(e)) for e in crosslinks})
            carbonyls = [co for (co, _n, _c) in backbone]
            offending = None
            for comp in comps:
                cs = [c for c in carbonyls if c in comp]
                if len(cs) > 1:
                    offending = (comp, cs)
                    break
            if offending is None:
                break
            comp, cs = offending
            path = _path(mol, cs[0], cs[1], comp)
            cut = None
            if path:
                for u, v in zip(path, path[1:]):
                    key = tuple(sorted((u, v)))
                    if key in bb_break or key in {tuple(sorted(e)) for e in crosslinks}:
                        continue
                    if mol.GetBondBetweenAtoms(u, v).IsInRing():
                        continue
                    cut = (u, v)
                    break
                if cut is None:
                    for u, v in zip(path, path[1:]):
                        cut = (u, v); break
            if cut is None:
                break
            crosslinks.append(cut)
        return crosslinks

    # ── matching ─────────────────────────────────────────────────────

    def match_one(self, capped_smiles):
        m = Chem.MolFromSmiles(capped_smiles)
        if m is None:
            return None, 0.0
        fp = AllChem.GetMorganFingerprintAsBitVect(m, self.radius,
                                                   nBits=self.nbits, useChirality=True)
        best_id, best = None, -1.0
        for bid, bfp in self._block_fp.items():
            s = DataStructs.TanimotoSimilarity(fp, bfp)
            if s > best:
                best, best_id = s, bid
        return best_id, best

    # ── public API ───────────────────────────────────────────────────

    def tokenize(self, smiles):
        """Full decomposition: SMILES -> ordered per-residue block matches."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Could not parse SMILES")
        backbone = self._backbone_amides(mol)
        crosslinks = self._detect_crosslinks(mol, backbone)
        frag_idx, frag_smi = self._capped_fragments(mol, backbone, crosslinks)

        # map atom -> fragment id
        atom_frag = {}
        for fi, idxs in enumerate(frag_idx):
            for a in idxs:
                atom_frag[a] = fi
        # residue adjacency: frag(co) -> frag(n) for each backbone amide
        nxt, has_pred = {}, set()
        carbonyl_frag = {}
        for (co, n, _c) in backbone:
            fc, fn = atom_frag.get(co), atom_frag.get(n)
            carbonyl_frag.setdefault(fc, co)
            if fc is not None and fn is not None and fc != fn:
                nxt[fc] = fn
                has_pred.add(fn)
        # order fragments by walking the chain
        order = []
        if nxt:
            start = next((f for f in nxt if f not in has_pred), next(iter(nxt)))
            cur, seen = start, set()
            while cur is not None and cur not in seen:
                order.append(cur); seen.add(cur)
                cur = nxt.get(cur)
            for fi in range(len(frag_idx)):
                if fi not in seen and fi in carbonyl_frag:
                    order.append(fi)
        else:
            order = [fi for fi in range(len(frag_idx))]

        residues = []
        n_debris = 0
        for fi in order:
            # skip pure cap/fragment with no backbone carbonyl
            if fi not in carbonyl_frag and len(frag_idx) > 1:
                continue
            capped = frag_smi[fi]
            # A crosslink cut can orphan a backbone carbonyl as a tiny capped
            # stub (e.g. CNC=O) that is not a real residue. The Ac/NHMe caps
            # alone are 5 heavy atoms, so anything below ~7 has no residue body.
            cm = Chem.MolFromSmiles(capped)
            if (len(frag_idx) > 1 and cm is not None
                    and cm.GetNumHeavyAtoms() < 7):
                n_debris += 1
                continue
            bid, sim = self.match_one(capped)
            residues.append({"capped": capped, "block": bid,
                             "tanimoto": round(sim, 3)})

        return {
            "n_residues": len(residues),
            "sequence": "-".join(r["block"] for r in residues if r["block"]),
            "residues": residues,
            "n_crosslinks": len(crosslinks),
            "n_debris": n_debris,
            "mean_tanimoto": round(
                sum(r["tanimoto"] for r in residues) / max(len(residues), 1), 3),
        }

    def roundtrip(self, sequence_str, assembler=None):
        """Assemble a ResToken sequence, decompose it, report block recovery."""
        from collections import Counter
        from restoken.src.cyclic_assembler import CyclicPeptideAssembler
        asm = assembler or CyclicPeptideAssembler()
        smi = asm.assemble(sequence_str)
        if smi is None:
            return {"ok": False, "reason": "assembly failed"}
        result = self.tokenize(smi)
        original = [s.strip() for s in sequence_str.split("-")]
        recovered = [r["block"] for r in result["residues"]]
        n_match = sum((Counter(recovered) & Counter(original)).values())
        return {
            "ok": True,
            "original": original,
            "recovered": recovered,
            "exact_multiset": Counter(recovered) == Counter(original),
            "recovery": round(n_match / len(original), 3),
            "mean_tanimoto": result["mean_tanimoto"],
            "n_residues": result["n_residues"],
        }


# ── graph helpers ──────────────────────────────────────────────────────

def _components(mol, broken):
    n = mol.GetNumAtoms()
    seen = [False] * n
    comps = []
    for s in range(n):
        if seen[s]:
            continue
        comp, dq = set(), deque([s]); seen[s] = True
        while dq:
            i = dq.popleft(); comp.add(i)
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if not seen[j] and tuple(sorted((i, j))) not in broken:
                    seen[j] = True; dq.append(j)
        comps.append(comp)
    return comps


def _path(mol, src, dst, allowed):
    prev = {src: None}
    dq = deque([src])
    while dq:
        cur = dq.popleft()
        if cur == dst:
            out = []
            while cur is not None:
                out.append(cur); cur = prev[cur]
            return out[::-1]
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j in allowed and j not in prev:
                prev[j] = cur; dq.append(j)
    return None
