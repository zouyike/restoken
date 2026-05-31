import glob, sys
import pyrosetta
from pyrosetta import init
from pyrosetta.rosetta.core.chemical import ChemicalManager
from pyrosetta.rosetta.core.pose import Pose

pdir = "/scratch/genesis/NCAA_tokenization/experiments/docking/params"
params = sorted(glob.glob(pdir + "/*.params"))

# Load all params together via -extra_res_fa
init("-beta_nov16 -mute all -extra_res_fa " + " ".join(params))

chm = ChemicalManager.get_instance()
rts = chm.residue_type_set("fa_standard")

ok, bad = [], []
for p in params:
    name = p.split("/")[-1].replace(".params", "")
    try:
        rt = rts.name_map(name)
        nat = rt.natoms()
        props = []
        if rt.is_polymer(): props.append("POLYMER")
        if rt.is_alpha_aa(): props.append("ALPHA_AA")
        if rt.is_beta_aa(): props.append("BETA_AA")
        lc = rt.lower_connect_id()
        uc = rt.upper_connect_id()
        ok.append((name, nat, lc, uc, " ".join(props)))
    except Exception as e:
        bad.append((name, str(e)[:80]))

print("=== LOADED OK ===")
for n, nat, lc, uc, pr in ok:
    print(f"{n}: natoms={nat} lower_conn={lc} upper_conn={uc} [{pr}]")
print("=== FAILED ===")
for n, e in bad:
    print(f"{n}: {e}")
print(f"\nSUMMARY: {len(ok)}/{len(params)} loaded")
