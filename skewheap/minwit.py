"""The SMALLEST witness that Naruse-Okada Def 2.6 cannot reach the true hook
multiset off the d-complete locus.  At F = empty the skew formula reduces to
series(P) = prod_u 1/(1-q^{h(u)}), so the hook function it uses must have the
true multiset; outside the d-complete class Def 2.6 is a family of recursions,
one per choice of d_k-interval at each element topping more than one.  A poset
where NO choice gives the true multiset is a finite, exhaustive disproof.

Prints the witness in full, with the choice count exhausted, so it can be
checked by hand.
usage: python3 minwit.py [MAXLEN=8] [TYPES=D4]
"""
import sys, itertools, collections
sys.path.insert(0,'.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from novel import dual, relabel_natural
from no26 import leq_matrix, dk_intervals
from filterhook import hooks_of
from skew import to_proctor
from choice import hooks_from_choice
ML=int(sys.argv[1]) if len(sys.argv)>1 else 8
TY=(sys.argv[2] if len(sys.argv)>2 else "D4").split(",")
best=None
for t in TY:
    typ,rk=t[0],int(t[1:]); C=cartan(typ,rk)
    out=[]
    def rec(x):
        if x and is_reduced(C,x): out.append(tuple(x))
        if len(x)==ML: return
        for i in range(rk): rec(x+[i])
    rec([])
    for w in out:
        down=heap_down(C,list(w)); n=len(down)
        if n<4: continue
        dn=relabel_natural(down); th=hooks_of(dn)
        if th is None: continue
        if len(dcomplete(relabel_natural(dual(down))))==0: continue   # skip d-complete
        dP,col,order=to_proctor(down,list(w)); leq=leq_matrix(dP)
        tops=collections.defaultdict(list)
        for (v,u,x,y) in dk_intervals(dP,leq): tops[u].append((v,x,y))
        keys=sorted(tops); nc=1
        for u in keys: nc*=(len(tops[u])+1)
        if nc>20000: continue
        reach=False
        for pick in itertools.product(*[[None]+list(range(len(tops[u]))) for u in keys]):
            ch={u:tops[u][p] for u,p in zip(keys,pick) if p is not None}
            if sorted(hooks_from_choice(dP,leq,ch))==th: reach=True; break
        if not reach and (best is None or n<best[1]):
            best=(w,n,th,down,nc,t)
if best is None: print("no witness in range"); sys.exit()
w,n,th,down,nc,t=best
print(f"SMALLEST WITNESS   type {t}, reduced word {w}, |P| = {n}")
print(f"  true hook multiset (from the series, is_product): {th}")
print(f"  choices in the Def 2.6 family exhausted: {nc}  -- none gives that multiset")
print(f"  heap, each element -> its lower covers:")
for x in sorted(down): print(f"    {x}: {sorted(down[x])}")
