"""Non-simply-laced heaps.  h(beta) = <beta^vee, mu> needs the bilinear form.

C_ij = <alpha_j, alpha_i^vee> = 2(alpha_i,alpha_j)/(alpha_i,alpha_i).
Set (alpha_i,alpha_i) = 2 eps_i; symmetry forces C_ij eps_i = C_ji eps_j.
Then (alpha_i,alpha_j) = C_ij eps_i, and for mu = sum m_i omega_i,
(alpha_i,mu) = m_i eps_i, so

   h(beta) = <beta^vee,mu> = 2 sum_i c_i eps_i m_i / (beta,beta).
"""
from fractions import Fraction as Fr

def cartan(typ,n):
    C=[[0]*n for _ in range(n)]
    for i in range(n): C[i][i]=2
    def link(i,j,a,b): C[i][j]=a; C[j][i]=b
    if typ=="A":
        for i in range(n-1): link(i,i+1,-1,-1)
    elif typ=="D":
        for i in range(n-2): link(i,i+1,-1,-1)
        link(n-3,n-1,-1,-1)
    elif typ=="B":   # alpha_{n-1} short
        for i in range(n-2): link(i,i+1,-1,-1)
        link(n-2,n-1,-2,-1)
    elif typ=="C":   # alpha_{n-1} long
        for i in range(n-2): link(i,i+1,-1,-1)
        link(n-2,n-1,-1,-2)
    elif typ=="G":   # G2
        C=[[2,-3],[-1,2]]
    elif typ=="F":   # F4
        C=[[2,-1,0,0],[-1,2,-2,0],[0,-1,2,-1],[0,0,-1,2]]
    else: raise ValueError(typ)
    return C

def eps(C):
    n=len(C); e=[None]*n; e[0]=Fr(1); seen=[0]
    while seen:
        i=seen.pop()
        for j in range(n):
            if C[i][j]!=0 and i!=j and e[j] is None:
                e[j]=e[i]*Fr(C[i][j],C[j][i]); seen.append(j)
    assert all(x is not None for x in e)
    for i in range(n):
        for j in range(n):
            assert C[i][j]*e[i]==C[j][i]*e[j], (i,j)
    return e

def adjacent(C,i,j): return i==j or C[i][j]!=0
def act(C,i,v):
    pair=sum(v[j]*C[i][j] for j in range(len(v)))
    w=list(v); w[i]-=pair; return tuple(w)
def roots_of_word(C,word):
    n=len(C); out=[]
    for k,ik in enumerate(word):
        v=[0]*n; v[ik]=1; v=tuple(v)
        for i in reversed(word[:k]): v=act(C,i,v)
        out.append(v)
    return out
def is_reduced(C,word):
    R=roots_of_word(C,word)
    return len(set(R))==len(R) and all(all(c>=0 for c in r) and sum(r)>0 for r in R)
def heap_down(C,word):
    L=len(word); d=[set() for _ in range(L)]
    for l in range(L):
        for k in range(l):
            if adjacent(C,word[k],word[l]): d[l]|=d[k]|{k}
    return [frozenset(s) for s in d]
def hooks(C,word,m):
    e=eps(C); out=[]
    for r in roots_of_word(C,word):
        bb=sum(Fr(r[i]*r[j])*C[i][j]*e[i] for i in range(len(C)) for j in range(len(C)))
        num=2*sum(Fr(r[i])*e[i]*m[i] for i in range(len(C)))
        h=num/bb
        out.append(h)
    return out
