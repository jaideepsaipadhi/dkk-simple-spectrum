"""The product-but-NOT-d-complete heaps: the class that is not explained by
Peterson-Proctor.  Dump the smallest ones with full structure, and test whether each is a
SLANT SUM (gluing two posets at a single element), since leaf posets are built from basic
leaf posets by slant/joint sums -- so a slant-decomposable example is probably already
covered by Ishikawa-Tagawa, while an indecomposable one is not.
usage: python3 novel.py [MAXLEN=9] [PMAX=14] [TYPES=D4]
"""
import sys, collections, itertools
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from crit import is_product
from wfast import Wvec_fast

def dual(down):
    n=len(down); up=[[] for _ in range(n)]
    for u in range(n):
        for v in down[u]: up[v].append(u)
    return [frozenset(up[u]) for u in range(n)]

def relabel_natural(down):
    n=len(down); need=[set(down[u]) for u in range(n)]
    order=[]; placed=set()
    while len(order)<n:
        for u in range(n):
            if u in placed: continue
            if need[u] <= placed: order.append(u); placed.add(u); break
        else: raise RuntimeError('cycle')
    pos={u:i for i,u in enumerate(order)}
    return [frozenset(pos[v] for v in down[u]) for u in order]

def ser(down, N=30):
    n=len(down); num=list(Wvec_fast(down))+[0]*(N+1); num=num[:N+1]
    for i in range(1,n+1):
        for d in range(i,N+1): num[d]+=num[d-i]
    return num[:N+1]

def closure(down):
    n=len(down); le=[[False]*n for _ in range(n)]
    for u in range(n): le[u][u]=True
    for _ in range(n):
        for u in range(n):
            for v in range(n):
                if le[u][v]:
                    for x in down[v]:
                        if not le[u][x]: le[u][x]=True
    return le   # le[u][v] True means v <= u

def covers(down):
    n=len(down); le=closure(down); cov=[[] for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v and le[u][v]:
                if not any(w!=u and w!=v and le[u][w] and le[w][v] for w in range(n)):
                    cov[u].append(v)
    return cov

def is_slant_sum(down):
    """does removing one element x disconnect the Hasse diagram into pieces that only
    meet at x?  (a cut vertex) -- the shape of a slant/joint sum."""
    n=len(down); cov=covers(down)
    adj={u:set() for u in range(n)}
    for u in range(n):
        for v in cov[u]: adj[u].add(v); adj[v].add(u)
    for x in range(n):
        rest=[u for u in range(n) if u!=x]
        if not rest: continue
        seen={rest[0]}; st=[rest[0]]
        while st:
            a=st.pop()
            for b in adj[a]:
                if b!=x and b not in seen: seen.add(b); st.append(b)
        if len(seen)!=len(rest): return True, x
    return False, None

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "D4").split(",")
    hits=[]
    for t in TY:
        typ, rk = t[0], int(t[1:]); C = cartan(typ, rk)
        out=[]
        def rec(w):
            if w and is_reduced(C, w): out.append(tuple(w))
            if len(w)==ML: return
            for i in range(rk): rec(w+[i])
        rec([])
        for w in out:
            down = heap_down(C, list(w)); n=len(down)
            if n<2 or n>PMAX: continue
            dn = relabel_natural(down)
            dc = (len(dcomplete(relabel_natural(dual(down))))==0)
            if dc: continue
            s = ser(dn)
            try: ok,bad,m = is_product(s,29,n)
            except AssertionError: continue
            if not (ok and sum(m.values())==n): continue
            hooks = sorted(k for k,v in m.items() for _ in range(v))
            sl,x = is_slant_sum(dn)
            hits.append((n,t,w,hooks,sl,x,dn,covers(dn)))
    hits.sort(key=lambda h:(h[0],h[1],len(h[2])))
    print(f"types {TY}, words to length {ML}, |P| <= {PMAX}")
    print(f"{len(hits)} heaps with product RPP that are NOT d-complete")
    nsl = sum(1 for h in hits if not h[4])
    print(f"   slant-decomposable (cut vertex): {len(hits)-nsl};  indecomposable: {nsl}\n")
    shown=0
    for n,t,w,hooks,sl,x,dn,cov in hits:
        if shown>=6: break
        print(f"|P|={n}  {t} w={''.join('s'+str(i) for i in w)}")
        print(f"    hooks = {hooks}   sum={sum(hooks)}  prod={1}")
        print(f"    cut vertex: {x if sl else 'NONE (indecomposable)'}")
        print(f"    covers (u: elements directly below u): "
              + "; ".join(f"{u}:{sorted(cov[u])}" for u in range(n) if cov[u]))
        shown+=1

if __name__ == "__main__":
    main()
