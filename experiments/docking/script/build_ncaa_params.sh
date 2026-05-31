#!/bin/bash
# Build a Rosetta POLYMER (.params) for one NCAA from its capped free-amino-acid SMILES.
# Capped form must be ACE-<residue>-NME, i.e. CH3-C(=O)-NH-<CA...>-C(=O)-NH-CH3,
# so write_M_session can detect the ACE/NME caps and the N-CA-C backbone.
#
# Validated chain (all on HPC, no MOE):
#   SMILES -> RDKit 3D -> openbabel mol2 -> write_M_session (py3, patched)
#   -> insert @<TRIPOS>SUBSTRUCTURE -> molfile_to_params_polymer_v4 (py2, --polymer)
#
# Usage: build_ncaa_params.sh <NAME> "<capped_SMILES>" [out_dir]
set -euo pipefail

NAME="$1"; SMILES="$2"
OUT="${3:-/scratch/genesis/NCAA_tokenization/experiments/docking/params}"

BIO=/public/home/genesis/miniconda3/envs/bio_env/bin
PY2=/public/home/genesis/miniconda3/envs/py2/bin/python2
RD=/public/home/genesis/RD/make_params
WMS=/scratch/genesis/NCAA_tokenization/experiments/docking/script/write_M_session_local.py
V4="$RD/molfile_to_params_polymer_v4.py"

export BABEL_LIBDIR=/public/home/genesis/miniconda3/envs/bio_env/lib/openbabel/3.1.0
export BABEL_DATADIR=/public/home/genesis/miniconda3/envs/bio_env/share/openbabel/3.1.0

mkdir -p "$OUT"; cd "$OUT"
# registry needed in cwd by v4 (auto-assigns TLC + appends); use a local copy
[ -f ncaa_registration.dat ] || cp "$RD/ncaa_registration.dat" ./ncaa_registration.dat

echo "[1/4] $NAME : SMILES -> 3D sdf (RDKit)"
$BIO/python - "$SMILES" "$NAME" <<'PY'
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
smi, name = sys.argv[1], sys.argv[2]
m = Chem.AddHs(Chem.MolFromSmiles(smi))
if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
    AllChem.EmbedMolecule(m, randomSeed=1, useRandomCoords=True)
AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
Chem.MolToMolFile(m, name + ".sdf")
print("   atoms=", m.GetNumAtoms())
PY

echo "[2/4] $NAME : sdf -> mol2 (openbabel)"
$BIO/obabel "$NAME.sdf" -O "$NAME.mol2" 2>/dev/null

echo "[3/4] $NAME : write_M_session + SUBSTRUCTURE marker"
$BIO/python "$WMS" "$NAME.mol2" --overwrite >/dev/null 2>&1 || true
/public/home/genesis/miniconda3/bin/python - "$NAME.mol2" <<'PY'
import sys
p = sys.argv[1]; L = open(p).read().splitlines(); out = []; done = False
for ln in L:
    if (not done) and ln.startswith('M  '):
        out.append('@<TRIPOS>SUBSTRUCTURE'); done = True
    out.append(ln)
open(p, 'w').write('\n'.join(out) + '\n')
assert done, "no M POLY records written -- backbone/caps not detected"
PY

echo "[4/4] $NAME : molfile_to_params_polymer_v4 --polymer (py2)"
rm -f "$NAME.params"
$PY2 "$V4" --overwrite --polymer --name "$NAME" "$NAME.mol2" 2>&1 | grep -iE "wrote|error|traceback|warning: very" || true

if grep -q "^TYPE POLYMER" "$NAME.params" 2>/dev/null; then
    echo "OK -> $OUT/$NAME.params  ($(grep -m1 '^PROPERTIES' "$NAME.params" | cut -c12-))"
else
    echo "FAIL: $NAME.params is not TYPE POLYMER"; exit 1
fi
