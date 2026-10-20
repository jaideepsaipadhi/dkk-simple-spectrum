"""The type-D criterion, transported: simple spectrum <=> the greedy top
truncation of the colour-t profile of I_kappa to total v_t is multiplicity
free.  In type D the profile is (1,2^{k-1},1) and this reduces to v_t <= 2.
In E6 the profile depends on kappa, so the threshold does too.
usage: python3 ecrit2.py [E6 D6 ...]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama

def prof_of(rs,k,t):
    pr=proj(rs,nakayama(rs)[k]); return {d:pr.dim.get((d,t),0) for d in range(pr.top+1)}
def greedy(prof,m):
    out={}
    for d in sorted(prof,reverse=True):
        if m<=0: break
        take=min(prof[d],m)
        if take: out[d]=take; m-=take
    return out,m
def threshold(prof):
    """largest m for which the greedy truncation is multiplicity free"""
    tot=sum(prof.values()); best=0
    for m in range(tot+1):
        g,left=greedy(prof,m)
        if left==0 and all(c<=1 for c in g.values()): best=m
        else: break
    return best

for arg in (sys.argv[1:] or ['D6','E6']):
    T,N=arg[0],int(arg[1:]); rs=RootSystem(T,N)
    t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    KR=[k for k in rs.I if 2<=k<=N-2] if T=='D' else list(rs.I)
    th={k:threshold(prof_of(rs,k,t)) for k in KR}
    print(f"\n{T}{N} (node {t}): thresholds v_t <= T(kappa):  {th}")
    tp=fp=fn=tn=0; ex=None
    for k in KR:
        lam=rs.fundamental(k)
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            _,rows,dv=soc_sigma(rs,k,w)
            simple=all(c<=1 for r in rows.values() for c in r.values())
            pred = dv[t-1]<=th[k]
            if pred and simple: tp+=1
            elif pred and not simple: fp+=1
            elif simple:
                fn+=1
                if ex is None: ex=(k,w,dv[t-1],th[k])
            else: tn+=1
    print(f"   rule 'v_t <= T(kappa)':  correct {tp+tn}/{tp+fp+fn+tn}"
          f"   false-positive {fp}   false-negative {fn}")
    if ex: print(f"   e.g. simple but above threshold: kappa={ex[0]} v_t={ex[2]} > T={ex[3]} w={ex[1]}")
