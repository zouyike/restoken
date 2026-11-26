import pandas as pd, numpy as np
WD="/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
df=pd.read_csv(f"{WD}/data/CycPeptMPDB_Peptide_All.csv", low_memory=False)
df['PAMPA']=pd.to_numeric(df['PAMPA'],errors='coerce')
d=df[df['PAMPA'].notna()][['CycPeptMPDB_ID','SMILES','PAMPA','Monomer_Length']].copy()
d=d.dropna(subset=['SMILES']).drop_duplicates(subset=['SMILES'])
# RDKit sanitize
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
def ok(s):
    m=Chem.MolFromSmiles(s); return m is not None
d=d[d['SMILES'].map(ok)].reset_index(drop=True)
rng=np.random.RandomState(42)
idx=rng.permutation(len(d)); ntest=int(len(d)*0.1)
te=d.iloc[idx[:ntest]].copy(); tr=d.iloc[idx[ntest:]].copy()
tr.to_csv(f"{WD}/data/pampa_train.csv",index=False)
te.to_csv(f"{WD}/data/pampa_test.csv",index=False)
print(f"valid={len(d)}  train={len(tr)}  test={len(te)}")
print(f"PAMPA mean/std (train): {tr.PAMPA.mean():.3f} / {tr.PAMPA.std():.3f}")
print(f"Monomer_Length: {d.Monomer_Length.min()}-{d.Monomer_Length.max()}, median {d.Monomer_Length.median()}")
