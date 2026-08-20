"""
Pi(Q) for Q Dynkin is self-injective, so P = Pi e_k is injective with simple
socle.  Conjecture: EVERY graded submodule N of P has N_t equal to a TOP
interval of (P)_t -- i.e. the quotient's colour-t multiplicities are always the
greedy BOTTOM truncation.  If so, Lemma T is a structural fact about P, not a
property of our particular quotient.

Enumerate ALL graded submodules and test.
"""
from itertools import product
from qvar import rref, P
from proj2 import Proj2
from dfs import subspaces
from rootsys import RootSystem

_S={}
def subs(m):
    if m not in _S: _S[m]=subspaces(m)
    return _S[m]

def all_submodules(rs,k):
    pr=Proj2(rs,k); top=pr.top
    out=[]
    def preimage(d,i,N):
        m=pr.dim.get((d,i),0)
        if m==0: return []
        eqs=[]
        for j in rs.I:
            nj=pr.dim.get((d+1,j),0)
            if nj==0 or (d,i,j) not in pr.arrow: continue
            R,piv=rref([r[:] for r in N.get((d+1,j),[])],nj)
            imgs=[]
            for b in range(m):
                e=[1 if x==b else 0 for x in range(m)]
                img=pr.act(d,i,e,j)
                for r,c in zip(R,piv):
                    if img[c]%P:
                        g=img[c]; img=[(img[x]-g*r[x])%P for x in range(nj)]
                imgs.append(img)
            for x in range(nj):
                row=[imgs[b][x]%P for b in range(m)]
                if any(y%P for y in row): eqs.append(row)
        if not eqs: return [[1 if x==b else 0 for x in range(m)] for b in range(m)]
        R,piv=rref(eqs,m); free=[c for c in range(m) if c not in piv]
        ker=[]
        for f in free:
            x=[0]*m; x[f]=1
            for r,c in zip(R,piv): x[c]=(-r[f])%P
            ker.append(x)
        return ker
    def rec(d,N):
        if d<0:
            out.append(dict(N)); return
        verts=[i for i in rs.I if pr.dim.get((d,i),0)]
        if not verts: rec(d-1,N); return
        pres={i:preimage(d,i,N) for i in verts}
        def walk(tt,N2):
            if tt==len(verts): rec(d-1,N2); return
            i=verts[tt]; basis=pres[i]; m=pr.dim[(d,i)]
            for sub in subs(len(basis)):
                rows=[[sum(co[u]*basis[u][c] for u in range(len(basis)))%P
                       for c in range(m)] for co in sub]
                R,_=rref(rows,m)
                N3=dict(N2); N3[(d,i)]=R
                walk(tt+1,N3)
        walk(0,N)
    rec(top,{})
    return pr,out

for typ,n,k in [("D",4,2),("D",5,2),("D",5,3),("D",6,2)]:
    rs=RootSystem(typ,n); t=n-2
    pr,subm=all_submodules(rs,k)
    degs=[d for d in range(pr.top+1) if pr.dim.get((d,t),0)]
    seq=[pr.dim[(d,t)] for d in degs]
    def bottom_trunc(total):
        o={}; 
        for s,d in zip(seq,degs):
            if total<=0: break
            take=min(s,total); o[d]=take; total-=take
        return o
    bad=0
    for N in subm:
        q={d: pr.dim[(d,t)]-len(N.get((d,t),[])) for d in degs}
        q={d:c for d,c in q.items() if c}
        if q != bottom_trunc(sum(q.values())): bad+=1
    print(f"D{n} e_{k}: {len(subm)} graded submodules; "
          f"{bad} whose colour-{t} quotient is NOT a bottom-truncation")
