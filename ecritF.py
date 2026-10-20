"""Candidate criterion for simple spectrum in type E, stated uniformly.

Let t be the trivalent node, arms the components of Gamma\{t}, and arm(kappa)
the arm containing kappa (empty if kappa=t).  Let j- and j+ be the two LONG
arms' far ends (E6: 1 and 6; E7: 1 and 7), and n- , n+ the arm neighbours of t
in those arms (3 and 5).  T(kappa) from the colour-t profile as before.  Then

  simple  <=>  v_t <= T(kappa)
               OR  ( T(kappa) = 2  and  v_t = 3
                     and v_{n-}, v_{n+} both even
                     and v_{j-} = [j- not in arm(kappa)] (mod 2)
                     and v_{j+} = [j+ not in arm(kappa)] (mod 2) )

Fitted on E6; E7 is out-of-sample.
usage: python3 ecritE.py E6|E7
"""
import sys
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
T_,N=sys.argv[1][0],int(sys.argv[1][1:])
rs=RootSystem(T_,N); t=4
# arms of Gamma minus t
arms=[]; seen=set()
for i in rs.I:
    if i==t or i in seen: continue
    st=[i]; c=[]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x); c.append(x); st+=[y for y in rs.adj[x] if y!=t and y not in seen]
    arms.append(sorted(c))
long_arms=[a for a in arms if len(a)>1]
def far(a):   # the end of arm a not adjacent to t
    return max(a, key=lambda x: min(len(a), abs(x)) if False else (0 if t in rs.adj[x] else 1)*100+x)
ends=[]; nbrs=[]
for a in long_arms:
    nb=[x for x in a if t in rs.adj[x]][0]
    fr=[x for x in a if t not in rs.adj[x]]
    ends.append(max(fr)); nbrs.append(nb)
ARM={}
for a in arms:
    for i in a: ARM[i]=set(a)
ARM[t]=set()
def Tof(k):
    pr=proj(rs,nakayama(rs)[k])
    p=[pr.dim[(d,t)] for d in sorted((d for (d,i) in pr.dim if i==t),reverse=True)
       if pr.dim[(d,t)]>0]
    a=0
    while a<len(p) and p[a]==1: a+=1
    return a+1 if a<len(p) and p[a]>=2 else a
TH={k:Tof(k) for k in rs.I}
print(f"{T_}{N}: node {t}, arms {arms}, long-arm neighbours {nbrs}, far ends {ends}")
print(f"   thresholds {TH}")
def band(k,v):
    if any(v[j-1]%2 for j in nbrs): return False
    for j in ends:
        if v[j-1]%2 != (0 if j in ARM[k] else 1): return False
    return True
def pred(k,v): return v[t-1]<=TH[k] or (v[t-1]==TH[k]+1 and band(k,v))
tp=fp=fn=tn=0; exfp=exfn=None
for k in rs.I:
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        s=all(c<=1 for r in rows.values() for c in r.values()); p=pred(k,dv)
        if p and s: tp+=1
        elif p and not s: fp+=1; exfp=exfp or (k,dv)
        elif s: fn+=1; exfn=exfn or (k,dv)
        else: tn+=1
n=tp+fp+fn+tn
print(f"   correct {tp+tn}/{n}   false-positive {fp}   false-negative {fn}")
if exfp: print(f"     FP e.g. kappa={exfp[0]} v={exfp[1]}")
if exfn: print(f"     FN e.g. kappa={exfn[0]} v={exfn[1]}")
