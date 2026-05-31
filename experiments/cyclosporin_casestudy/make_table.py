from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import csv, os
# free amino acid SMILES at the Abu-2 position
rows = [
 ("Abu-2 (native)","ethyl","CCC(N)C(=O)O","0.00","0.00","native ref"),
 ("allylglycine","allyl (–CH2CH=CH2)","C=CCC(N)C(=O)O","-13.52","+0.56","MOE++, Rosetta neutral"),
 ("methionine","–CH2CH2SCH3","CSCCC(N)C(=O)O","-4.87","+0.45","MOE+, Rosetta neutral"),
 ("norleucine","n-butyl","CCCCC(N)C(=O)O","-3.53","-0.21","both engines improve"),
]
abu_mw=None
out=[]
for name,sc,smi,deint,ddg,note in rows:
    m=Chem.MolFromSmiles(smi)
    mw=Descriptors.MolWt(m); logp=Crippen.MolLogP(m)
    if abu_mw is None: abu_mw=mw; abu_logp=logp
    out.append(dict(analogue=name, side_chain=sc, freeAA_SMILES=smi,
                    freeAA_MW=round(mw,1), dMW_vs_Abu=round(mw-abu_mw,1),
                    cLogP=round(logp,2), dcLogP=round(logp-abu_logp,2),
                    dEint_MOE_kcal=deint, ddG_Rosetta_REU=ddg, verdict=note))
D="/scratch/genesis/NCAA_tokenization/experiments/cyclosporin_casestudy"
csvp=os.path.join(D,"csa_abu2_analogue_table.csv")
with open(csvp,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0].keys())); w.writeheader()
    for r in out: w.writerow(r)
print("wrote",csvp,"\n")
# pretty print
hdr=["analogue","side_chain","freeAA_MW","dMW","cLogP","dcLogP","dEint(MOE)","ddG(Ros)","verdict"]
print("  ".join(f"{h:>12}" for h in hdr))
for r in out:
    print(f"{r['analogue']:>12.12}  {r['side_chain']:>12.12}  {r['freeAA_MW']:>12}  {r['dMW_vs_Abu']:>12}  {r['cLogP']:>12}  {r['dcLogP']:>12}  {r['dEint_MOE_kcal']:>12}  {r['ddG_Rosetta_REU']:>12}  {r['verdict']}")
