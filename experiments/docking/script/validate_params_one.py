import sys
import pyrosetta
from pyrosetta import init
from pyrosetta.rosetta.core.chemical import ChemicalManager

p = sys.argv[1]
name = p.split("/")[-1].replace(".params", "")
try:
    init("-beta_nov16 -mute all -extra_res_fa " + p)
    chm = ChemicalManager.get_instance()
    rts = chm.residue_type_set("fa_standard")
    rt = rts.name_map(name)
    props = []
    if rt.is_polymer(): props.append("POLYMER")
    if rt.is_alpha_aa(): props.append("ALPHA_AA")
    if rt.is_beta_aa(): props.append("BETA_AA")
    print(f"OK {name}: natoms={rt.natoms()} lc={rt.lower_connect_id()} uc={rt.upper_connect_id()} [{' '.join(props)}]")
except Exception as e:
    print(f"FAIL {name}: {str(e)[:120]}")
