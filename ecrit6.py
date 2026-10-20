"""A candidate COMPLETE criterion for simple spectrum in E6.

  simple  <=>  v_t <= T(kappa)
               OR  ( v_t = T(kappa)+1  and  the band condition )

band condition, for the trivalent node t=4 of E6 with arms
A1={1,3}, A2={2}, A3={5,6} and arm(kappa) the one containing kappa (empty if
kappa=t):
   v_3 and v_5 both even, and
   v_1 = [1 not in arm(kappa)],  v_6 = [6 not in arm(kappa)]   (mod 2).
usage: python3 ecrit6.py
"""
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
rs=RootSystem('E',6); t=4
ARM={1:{1,3},3:{1,3},2:{2},5:{5,6},6:{5,6},4:set()}
def Tof(k):
    pr=proj(rs,nakayama(rs)[k])
    p=[pr.dim[(d,t)] for d in sorted((d for (d,i) in pr.dim if i==t),reverse=True)
       if pr.dim[(d,t)]>0]
    a=0
    while a<len(p) and p[a]==1: a+=1
    return a+1 if a<len(p) and p[a]>=2 else a
TH={k:Tof(k) for k in rs.I}
def band(k,v):
    if v[2]%2 or v[4]%2: return False
    return v[0]%2==(0 if 1 in ARM[k] else 1) and v[5]%2==(0 if 6 in ARM[k] else 1)
def pred(k,v): return v[t-1]<=TH[k] or (v[t-1]==TH[k]+1 and band(k,v))
tp=fp=fn=tn=0; exfp=exfn=None
for k in rs.I:
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        s=all(c<=1 for r in rows.values() for c in r.values()); p=pred(k,dv)
        if p and s: tp+=1
        elif p and not s:
            fp+=1; exfp=exfp or (k,dv,w)
        elif s:
            fn+=1; exfn=exfn or (k,dv,w)
        else: tn+=1
n=tp+fp+fn+tn
print(f"E6 thresholds T(kappa) = {TH}")
print(f"criterion correct {tp+tn}/{n}   false-positive {fp}   false-negative {fn}")
if exfp: print(f"  FP e.g. kappa={exfp[0]} v={exfp[1]} w={exfp[2]}")
if exfn: print(f"  FN e.g. kappa={exfn[0]} v={exfn[1]} w={exfn[2]}")
