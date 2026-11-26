import os, sys, json, argparse, numpy as np, pandas as pd
import torch, torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from transformers import RoFormerModel
from scipy.stats import spearmanr, pearsonr

WD = "/scratch/genesis/NCAA_tokenization/experiments/permeability_validation"
MODEL_DIR = f"{WD}/model"
sys.path.insert(0, MODEL_DIR)
from tokenizer.my_tokenizers import SMILES_SPE_Tokenizer

def get_tok():
    return SMILES_SPE_Tokenizer(f"{MODEL_DIR}/tokenizer/new_vocab.txt",
                                f"{MODEL_DIR}/tokenizer/new_splits.txt")

class DS(Dataset):
    def __init__(self, smiles, y, tok, maxlen):
        self.s=list(smiles); self.y=None if y is None else np.asarray(y,dtype=np.float32)
        self.tok=tok; self.maxlen=maxlen
    def __len__(self): return len(self.s)
    def __getitem__(self, i):
        e=self.tok(self.s[i], truncation=True, max_length=self.maxlen)
        item={"input_ids":torch.tensor(e["input_ids"]),
              "attention_mask":torch.tensor(e["attention_mask"])}
        if self.y is not None: item["y"]=torch.tensor(self.y[i])
        return item

def collate(batch, pad_id):
    m=max(len(b["input_ids"]) for b in batch)
    ids=torch.full((len(batch),m),pad_id,dtype=torch.long)
    am=torch.zeros((len(batch),m),dtype=torch.long)
    ys=[]
    for j,b in enumerate(batch):
        L=len(b["input_ids"]); ids[j,:L]=b["input_ids"]; am[j,:L]=b["attention_mask"]
        if "y" in b: ys.append(b["y"])
    out={"input_ids":ids,"attention_mask":am}
    if ys: out["y"]=torch.stack(ys)
    return out

class Regressor(nn.Module):
    def __init__(self, mdir, drop=0.1):
        super().__init__()
        self.enc=RoFormerModel.from_pretrained(mdir)
        h=self.enc.config.hidden_size
        self.head=nn.Sequential(nn.Dropout(drop), nn.Linear(h,256), nn.GELU(),
                                nn.Dropout(drop), nn.Linear(256,1))
    def forward(self, input_ids, attention_mask):
        o=self.enc(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
        mask=attention_mask.unsqueeze(-1).float()
        pooled=(o*mask).sum(1)/mask.sum(1).clamp(min=1e-6)   # mean-pool
        return self.head(pooled).squeeze(-1)

def evaluate(model, dl, dev, mu, sd):
    model.eval(); P=[]; Y=[]
    with torch.no_grad():
        for b in dl:
            p=model(b["input_ids"].to(dev), b["attention_mask"].to(dev)).cpu().numpy()
            P.append(p); Y.append(b["y"].numpy())
    P=np.concatenate(P)*sd+mu; Y=np.concatenate(Y)*sd+mu
    ss_res=((Y-P)**2).sum(); ss_tot=((Y-Y.mean())**2).sum()
    return {"r2":1-ss_res/ss_tot, "rmse":float(np.sqrt(((Y-P)**2).mean())),
            "spearman":float(spearmanr(Y,P).correlation),
            "pearson":float(pearsonr(Y,P)[0])}, P, Y

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--epochs",type=int,default=40)
    ap.add_argument("--bs",type=int,default=32)
    ap.add_argument("--lr",type=float,default=2e-5)
    ap.add_argument("--maxlen",type=int,default=512)
    ap.add_argument("--patience",type=int,default=6)
    a=ap.parse_args()
    dev="cuda" if torch.cuda.is_available() else "cpu"
    print("device:",dev, torch.cuda.get_device_name(0) if dev=="cuda" else "")
    tok=get_tok(); pad=tok.pad_token_id
    tr=pd.read_csv(f"{WD}/data/pampa_train.csv"); te=pd.read_csv(f"{WD}/data/pampa_test.csv")
    rng=np.random.RandomState(0); idx=rng.permutation(len(tr)); nv=int(len(tr)*0.1)
    va=tr.iloc[idx[:nv]]; trn=tr.iloc[idx[nv:]]
    mu=trn.PAMPA.mean(); sd=trn.PAMPA.std()
    print(f"train={len(trn)} val={len(va)} test={len(te)} | mu={mu:.3f} sd={sd:.3f}")
    mk=lambda d: DataLoader(DS(d.SMILES,(d.PAMPA-mu)/sd,tok,a.maxlen), batch_size=a.bs,
                            shuffle=False, collate_fn=lambda b:collate(b,pad))
    dl_tr=DataLoader(DS(trn.SMILES,(trn.PAMPA-mu)/sd,tok,a.maxlen),batch_size=a.bs,
                     shuffle=True, collate_fn=lambda b:collate(b,pad))
    dl_va=mk(va); dl_te=mk(te)
    model=Regressor(MODEL_DIR).to(dev)
    opt=torch.optim.AdamW(model.parameters(),lr=a.lr,weight_decay=0.01)
    lossfn=nn.MSELoss()
    best=-1e9; bad=0
    for ep in range(a.epochs):
        model.train(); tl=0
        for b in dl_tr:
            opt.zero_grad()
            p=model(b["input_ids"].to(dev),b["attention_mask"].to(dev))
            loss=lossfn(p,b["y"].to(dev)); loss.backward(); opt.step(); tl+=loss.item()
        vm,_,_=evaluate(model,dl_va,dev,mu,sd)
        print(f"ep{ep+1} trloss={tl/len(dl_tr):.4f} val_r2={vm['r2']:.3f} val_rho={vm['spearman']:.3f}",flush=True)
        if vm["spearman"]>best:
            best=vm["spearman"]; bad=0
            torch.save({"model":model.state_dict(),"mu":mu,"sd":sd,"maxlen":a.maxlen}, f"{WD}/output/perm_model.pt")
        else:
            bad+=1
            if bad>=a.patience: print("early stop"); break
    ck=torch.load(f"{WD}/output/perm_model.pt",map_location=dev); model.load_state_dict(ck["model"])
    tm,P,Y=evaluate(model,dl_te,dev,mu,sd)
    print("TEST:",json.dumps(tm,indent=2))
    json.dump(tm,open(f"{WD}/output/test_metrics.json","w"),indent=2)
    pd.DataFrame({"smiles":te.SMILES,"y_true":Y,"y_pred":P}).to_csv(f"{WD}/output/test_predictions.csv",index=False)

if __name__=="__main__": main()
