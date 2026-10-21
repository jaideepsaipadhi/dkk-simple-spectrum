"""Proctor's d-completeness, implemented from the definition:
 D_{m,n} (m,n>=2) = diamond + chain of m-2 above its top + chain of n-2 below its bottom, m+n elements.
 d_k-interval  [w,z] = D_{k-1,k-1};   d_k^- -set [w,..] = D_{k-2,k-1}  (d_3^- = diamond minus top).
 d_k-complete axioms:
  (A1) S a d_k^- set => exists z with [w,z] a d_k-interval
  (A2) [w,z] a d_k-interval => z covers nothing outside [w,z]
  (A3) [w,z] a d_k-interval => no w' != w with [w',z] a d_k-interval
 d-complete = d_k-complete for all k >= 3.
"""
import sys, itertools; sys.path.insert(0,'.')
def rel(down):
    n=len(down)
    le=[[False]*n for _ in range(n)]
    for u in range(n):
        le[u][u]=True
        for v in down[u]: le[v][u]=True
    return le
def covers(down):
    n=len(down); uc=[[] for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u in down[v] and not any(u in down[t] and t in down[v] for t in range(n)):
                uc[u].append(v)
    return uc
def shape(down, S):
    """If the induced poset on S is some D_{m,n} with a unique minimum, return (m,n,w,pair,top_or_None)."""
    le=rel(down); S=sorted(S)
    inc=[(a,b) for i,a in enumerate(S) for b in S[i+1:] if not le[a][b] and not le[b][a]]
    if len(inc)!=1: return None
    a1,a2=inc[0]
    L=[u for u in S if u!=a1 and u!=a2 and le[u][a1] and le[u][a2]]
    U=[u for u in S if u!=a1 and u!=a2 and le[a1][u] and le[a2][u]]
    if len(L)+len(U)+2!=len(S): return None
    def ischain(X): return all(le[x][y] or le[y][x] for x in X for y in X)
    if not (ischain(L) and ischain(U)): return None
    if not L: return None
    mins=[u for u in S if not any(v!=u and le[v][u] for v in S)]
    if len(mins)!=1: return None
    n_=len(L)+1; m_=len(U)+1          # |L| = n-1, |U| = m-1
    top=None
    if U:
        tops=[u for u in U if all(le[v][u] for v in U)]
        top=tops[0] if len(tops)==1 else None
    return (m_,n_,mins[0],(a1,a2),top)
def convex_sets(down, size):
    """all convex subsets of given size (small posets only)"""
    n=len(down); le=rel(down); out=[]
    for S in itertools.combinations(range(n), size):
        Ss=set(S); ok=True
        for u in S:
            for v in S:
                if le[u][v]:
                    for t in range(n):
                        if le[u][t] and le[t][v] and t not in Ss: ok=False; break
                if not ok: break
            if not ok: break
        if ok: out.append(S)
    return out
def interval(down, w, z):
    le=rel(down)
    if not le[w][z]: return None
    return [u for u in range(len(down)) if le[w][u] and le[u][z]]
def dcomplete(down, KMAX=None):
    n=len(down); le=rel(down); uc=covers(down)
    KMAX = KMAX if KMAX else (n+3)//2+1
    fails=[]
    for k in range(3, KMAX+1):
        szm = (k-2)+(k-3)+2          # |D_{k-2,k-1}| = (k-2)+(k-1) = 2k-3
        szf = (k-1)+(k-1)            # |D_{k-1,k-1}| = 2k-2
        # (A1)
        for S in convex_sets(down, szm):
            sh=shape(down,S)
            if sh is None: continue
            m_,n_,w,pair,top = sh
            if (m_,n_)!=(k-2,k-1): continue
            okz=False
            for z in range(n):
                I=interval(down,w,z)
                if I is None or len(I)!=szf: continue
                sh2=shape(down,I)
                if sh2 and (sh2[0],sh2[1])==(k-1,k-1) and set(S)<=set(I): okz=True; break
            if not okz: fails.append(("A1",k,S))
        # (A2),(A3)
        ints=[]
        for w in range(n):
            for z in range(n):
                I=interval(down,w,z)
                if I is None or len(I)!=szf: continue
                sh2=shape(down,I)
                if sh2 and (sh2[0],sh2[1])==(k-1,k-1) and sh2[2]==w: ints.append((w,z,I))
        for w,z,I in ints:
            if any(c not in I for c in uc[z] if False): pass
            below=[v for v in range(n) if z in uc[v]]      # elements z covers
            if any(v not in I for v in below): fails.append(("A2",k,w,z))
        byz={}
        for w,z,I in ints: byz.setdefault(z,[]).append(w)
        for z,ws in byz.items():
            if len(set(ws))>1: fails.append(("A3",k,z,sorted(set(ws))))
    return fails
if __name__ == "__main__":
    tests={
      "chain of 3":      [frozenset(),frozenset({0}),frozenset({0,1})],
      "V (w<a,b)":       [frozenset(),frozenset({0}),frozenset({0})],
      "Lambda (a,b<z)":  [frozenset(),frozenset(),frozenset({0,1})],
      "diamond":         [frozenset(),frozenset({0}),frozenset({0}),frozenset({0,1,2})],
      "diamond+tail up": [frozenset(),frozenset({0}),frozenset({0}),frozenset({0,1,2}),frozenset({0,1,2,3})],
    }
    print("VALIDATION (expect: chain pass, V fail A1, diamond pass):")
    for name,d in tests.items():
        f=dcomplete(d)
        print(f"  {name:18} -> {'d-complete' if not f else 'NOT d-complete: '+str(f[:2])}")
