import subprocess, re, sys, statistics, os
PY_BIN="/public/home/genesis/miniconda3/envs/SE3nv/bin/python"
JOBS=[("parent","input/parent_complex.pdb",[]),
      ("A86","input/analogue_p4_A86_complex.pdb",["--params","params_analogue/A86.params"]),
      ("A19","input/analogue_p9_A19_complex.pdb",["--params","params_analogue/A19.params"])]
NREP=3
out=open("rosetta/ddg_results.csv","w"); out.write("label,rep,dG,dSASA,closure\n"); out.flush()
agg={}
for label,cx,extra in JOBS:
    vals=[]
    for r in range(NREP):
        cmd=[PY_BIN,"script/interface_ddg.py","--complex",cx,"--params-dir","params_named",
             "--cyclize","--relax","--relax-rounds","1"]+extra
        p=subprocess.run(cmd,capture_output=True,text=True)
        t=p.stdout
        dG=re.search(r"dG_separated\s*:\s*(-?[\d.]+)",t)
        ds=re.search(r"dSASA_interface\s*:\s*([\d.]+)",t)
        cl=re.search(r"closure_after_fr\s*:\s*([\d.]+)",t)
        dG=float(dG.group(1)) if dG else float("nan")
        ds=float(ds.group(1)) if ds else float("nan")
        cl=float(cl.group(1)) if cl else float("nan")
        out.write(f"{label},{r},{dG:.3f},{ds:.1f},{cl:.3f}\n"); out.flush()
        vals.append(dG)
        print(f"{label} rep{r}: dG={dG:.2f} closure={cl:.3f}",flush=True)
    agg[label]=vals
out.close()
# summary + ddG vs parent
pm=statistics.mean(agg["parent"])
print("\n=== SUMMARY (mean dG, REU) ===")
for label in agg:
    m=statistics.mean(agg[label]); s=statistics.pstdev(agg[label])
    print(f"{label:8s} dG={m:8.2f} +/- {s:4.2f}   ddG_vs_parent={m-pm:+7.2f}")
