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

# Acetyl cap SMARTS (CH3-C(=O)-N): the match's 4th atom is the backbone N,
# used to locate where to transplant an N-modification onto a capped block.
_AC_SMARTS = Chem.MolFromSmarts("[CH3]C(=O)[N]")


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
                    # The amide's own carbonyl is a hard barrier: walking back
                    # into it reaches the *previous* residue's atoms and can find
                    # a backbone carbonyl on the wrong side (e.g. a pendant
                    # dimethylamide side chain whose N would otherwise tunnel
                    # through its carbonyl to the real backbone). Skip it.
                    if j == exclude_co:
                        continue
                    dist[j] = 1
                    if j in carbonyls:
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

    # ── N-modification recognition ───────────────────────────────────

    def _classify_nmod(self, mol, n_idx, co_prev, co_next):
        """Classify a residue's backbone nitrogen modification.

        The backbone N carries one of: a hydrogen (unmodified, an H-bond
        donor), an N-alkyl group (N-methyl/N-ethyl — the permeability-driving
        modification of drugs like LUNA18, which removes the donor), or it is
        part of a ring (proline-like, cyclic). N-alkylation is the single most
        important backbone feature a tokenizer must preserve, yet it is encoded
        by *which block* is chosen — and the library carries only a handful of
        N-alkyl variants. So we detect it structurally and report it per
        residue, independent of the matched block.

        Returns (kind, k): kind in {"NH","NME","NALK","NCY"}; k = number of
        carbons in the N-alkyl chain (1 = methyl, 2 = ethyl, ...).
        """
        a = mol.GetAtomWithIdx(n_idx)
        # Every backbone N sits in the head-to-tail macrocycle, so plain ring
        # membership is uninformative. Proline-like (cyclic) N is distinguished
        # by also belonging to a *small* ring (its pyrrolidine, size <= 8).
        if _in_small_ring(mol, n_idx, max_size=8):
            return ("NCY", 0)
        if a.GetTotalNumHs() > 0:
            return ("NH", 0)
        # Tertiary N: neighbours are the previous carbonyl, this residue's
        # alpha carbon, and the N-alkyl modifier. Identify the modifier as the
        # carbon neighbour that is neither the previous carbonyl nor the alpha
        # (the alpha is the one from which this residue's own carbonyl is
        # reachable without passing back through N).
        carbons = [nb.GetIdx() for nb in a.GetNeighbors()
                   if nb.GetAtomicNum() == 6 and nb.GetIdx() != co_prev]
        alpha = None
        for c in carbons:
            if _reaches(mol, c, co_next, blocked={n_idx}):
                alpha = c
                break
        modifier = [c for c in carbons if c != alpha]
        if not modifier:
            return ("NH", 0)
        k = _branch_carbon_count(mol, modifier[0], blocked={n_idx})
        return ("NME" if k == 1 else "NALK", k)

    def _set_backbone_nmod(self, capped_smiles, kind, k):
        """Return a capped residue whose backbone N carries the given
        modification, transplanting it onto a matched block's scaffold.

        ``capped_smiles`` is a library block's Ac-/NHMe-capped aa_smiles. We
        rewrite only its backbone nitrogen (the acetyl-capped N) so the side
        chain / scaffold stay the block's, while the N-modification is the one
        the real residue actually had. NCY (proline) residues are left intact —
        only proline blocks match them, and the ring is intrinsic.
        """
        if kind == "NCY":
            return capped_smiles
        m = Chem.MolFromSmiles(capped_smiles)
        if m is None:
            return capped_smiles
        ac = m.GetSubstructMatch(_AC_SMARTS)
        if not ac:
            return capped_smiles
        ac_co, n_bb = ac[1], ac[3]
        if m.GetAtomWithIdx(n_bb).IsInRing():
            return capped_smiles  # proline-like scaffold; leave the ring N
        rw = Chem.RWMol(m)
        # Strip any existing N-alkyl modifier on the backbone N: a carbon
        # neighbour (other than the acetyl carbonyl) that is a short terminal
        # branch carrying no carbonyl. The alpha carbon leads into the residue
        # body (and on to the C-terminal cap), so it is excluded by both tests.
        nbrs = [nb.GetIdx() for nb in rw.GetAtomWithIdx(n_bb).GetNeighbors()]
        modifier_atoms = []
        for c in nbrs:
            if c == ac_co or rw.GetAtomWithIdx(c).GetAtomicNum() != 6:
                continue
            branch = _branch_atoms(rw, c, blocked={n_bb})
            # alpha branch is large (carries the C-terminal cap); modifier is a
            # short alkyl with no carbonyl
            has_co = any(self._is_carbonyl_c(rw, b) for b in branch)
            if not has_co and len(branch) <= 3:
                modifier_atoms.append((c, branch))
        for _c, branch in modifier_atoms:
            for idx in sorted(branch, reverse=True):
                rw.RemoveAtom(idx)
        # re-resolve backbone N index after removals
        m2 = rw.GetMol()
        try:
            Chem.SanitizeMol(m2)
        except Exception:
            return capped_smiles
        rw = Chem.RWMol(m2)
        ac = rw.GetMol().GetSubstructMatch(_AC_SMARTS)
        if not ac:
            return capped_smiles
        n_bb = ac[3]
        if kind == "NH":
            return Chem.MolToSmiles(rw.GetMol())
        # attach a linear alkyl of k carbons
        prev = n_bb
        for _ in range(max(k, 1)):
            c = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(prev, c, Chem.BondType.SINGLE)
            prev = c
        out = rw.GetMol()
        try:
            Chem.SanitizeMol(out)
        except Exception:
            return capped_smiles
        return Chem.MolToSmiles(out)

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
        # residue backbone N: the residue owning carbonyl co2 has incoming
        # amino nitrogen n (bonded to the previous carbonyl co). Record both so
        # we can read off each residue's N-modification from the intact molecule.
        frag_backbone_n = {}
        for (co, n, co2) in backbone:
            fc, fn = atom_frag.get(co), atom_frag.get(n)
            carbonyl_frag.setdefault(fc, co)
            fco2 = atom_frag.get(co2)
            if fco2 is not None:
                frag_backbone_n[fco2] = (n, co, co2)
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
            nmod = ("NH", 0)
            if fi in frag_backbone_n:
                n_idx, co_prev, co_next = frag_backbone_n[fi]
                nmod = self._classify_nmod(mol, n_idx, co_prev, co_next)
            kind, k = nmod
            nmod_str = kind if kind != "NALK" else f"N-C{k}"
            residues.append({"capped": capped, "block": bid,
                             "tanimoto": round(sim, 3),
                             "nmod": kind, "nmod_k": k, "nmod_label": nmod_str})

        n_nmod = sum(1 for r in residues if r["nmod"] in ("NME", "NALK", "NCY"))
        return {
            "n_residues": len(residues),
            "sequence": "-".join(r["block"] for r in residues if r["block"]),
            "residues": residues,
            "n_crosslinks": len(crosslinks),
            "n_debris": n_debris,
            "n_nmod": n_nmod,
            "mean_tanimoto": round(
                sum(r["tanimoto"] for r in residues) / max(len(residues), 1), 3),
        }

    def reconstruct(self, result, assembler=None, nmod_correct=True):
        """Reassemble a tokenized result into a cyclic-peptide SMILES.

        With ``nmod_correct`` (default), each residue is built from its matched
        block's *scaffold* but its backbone nitrogen is set to the
        N-modification the real residue actually carried (detected in
        ``tokenize``). This preserves N-methylation / N-alkylation — LUNA18's
        permeability-defining feature — even though the library lacks a matching
        N-alkyl block, so the reconstruction's H-bond-donor count is not
        inflated by spurious backbone N-H. Without correction, the plain block
        round-trip is used.

        Returns the reconstructed SMILES, or None on failure.
        """
        from restoken.src.cyclic_assembler import CyclicPeptideAssembler
        asm = assembler or CyclicPeptideAssembler()
        residues = [r for r in result["residues"] if r.get("block")]
        if not residues:
            return None
        if not nmod_correct:
            return asm.assemble("-".join(r["block"] for r in residues))
        capped = []
        for r in residues:
            block_aa = asm._blocks[r["block"]]["structure"]["aa_smiles"]
            capped.append(self._set_backbone_nmod(
                block_aa, r["nmod"], r["nmod_k"]))
        return asm.assemble_advanced(capped, crosslinks=None)

    def crosslink_reencodability(self, smiles):
        """Report whether a molecule's side-chain bridges can be re-encoded by the
        assembler's crosslink primitive (assemble_advanced).

        The primitive re-forms a bridge as a single bond between two *side-chain*
        atoms of otherwise head-to-tail residues. That model holds only when each
        detected crosslink is atom-disjoint from the backbone amides: a bridge
        endpoint that is itself a backbone-amide atom means the macrocycle and the
        bridge share atoms (a fused polycycle, e.g. MK-0616), which cannot be
        expressed as head-to-tail residues + a disjoint staple.

        Returns dict: n_crosslinks, n_disjoint (re-encodable), n_fused (sharing a
        backbone atom), per-crosslink detail, and a boolean ``reencodable``.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Could not parse SMILES")
        backbone = self._backbone_amides(mol)
        crosslinks = self._detect_crosslinks(mol, backbone)
        bb_atoms = set()
        for co, n, _c in backbone:
            bb_atoms.add(co)
            bb_atoms.add(n)
        detail, n_fused = [], 0
        for u, v in crosslinks:
            shared = [a for a in (u, v) if a in bb_atoms]
            fused = bool(shared)
            n_fused += int(fused)
            detail.append({"bond": (u, v), "fused": fused, "shared_atoms": shared})
        n_disjoint = len(crosslinks) - n_fused
        return {
            "n_crosslinks": len(crosslinks),
            "n_disjoint": n_disjoint,
            "n_fused": n_fused,
            "detail": detail,
            "reencodable": len(crosslinks) > 0 and n_fused == 0,
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


def _in_small_ring(mol, idx, max_size=8):
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) <= max_size and idx in ring:
            return True
    return False


def _reaches(mol, start, target, blocked):
    """True if target is reachable from start without crossing a blocked atom."""
    if start == target:
        return True
    seen = set(blocked); seen.add(start)
    dq = deque([start])
    while dq:
        cur = dq.popleft()
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j == target:
                return True
            if j not in seen:
                seen.add(j); dq.append(j)
    return False


def _branch_atoms(mol, start, blocked):
    """All atom indices in the branch rooted at start, not crossing blocked."""
    seen = set(blocked)
    comp, dq = set(), deque([start])
    seen.add(start)
    while dq:
        cur = dq.popleft(); comp.add(cur)
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j not in seen:
                seen.add(j); dq.append(j)
    return comp


def _branch_carbon_count(mol, start, blocked):
    """Number of carbon atoms in the branch rooted at start."""
    return sum(1 for i in _branch_atoms(mol, start, blocked)
               if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)


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
