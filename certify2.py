"""RPP series of the witness heap in BOTH orientations, against the claimed
hook product.  DKK Def 5.2 has a < b when a > b as word positions, so later
letters are LOWER; the opposite choice gives the order dual, whose RPP series
differs.  Computed from scratch.
usage: python3 certify2.py
"""
D=26
W=[1,0,1,2,3,1,0,2]; ADJ={0:{1},1:{0,2,3},2:{1},3:{1}}
n=len(W)
def closure(rel):
    r=[row[:] for row in rel]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if r[i][k] and r[k][j]: r[i][j]=True
    return r
# orientation A: earlier position lower
A=[[False]*n for _ in range(n)]
for b in range(n):
    for a in range(b):
        if W[a]==W[b] or W[b] in ADJ[W[a]]: A[a][b]=True
A=closure(A)
# orientation B (DKK Def 5.2): later position lower
B=[[False]*n for _ in range(n)]
for b in range(n):
    for a in range(b):
        if W[a]==W[b] or W[b] in ADJ[W[a]]: B[b][a]=True
B=closure(B)
def mul(a,b):
    c=[0]*(D+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y and i+j<=D: c[i+j]+=x*y
    return c
def geom(k):
    z=[0]*(D+1)
    for i in range(0,D+1,k): z[i]=1
    return z
def rpp_series(less):
    filt=[]
    for m in range(1<<n):
        S={i for i in range(n) if m>>i&1}
        if all(not(less[x][y] and x in S and y not in S) for x in range(n) for y in range(n)):
            filt.append(frozenset(S))
    ne=[F for F in filt if F]
    g={}
    for F in sorted(ne,key=len):
        acc=[0]*(D+1); acc[0]=1
        for G in ne:
            if G<F:
                for i in range(D+1): acc[i]+=g[G][i]
        k=len(F); sh=[0]*(D+1); sh[k]=1        # q^{|F|}, one coefficient
        g[F]=mul(mul(sh,geom(k)),acc)
    s=[0]*(D+1); s[0]=1
    for F in ne:
        for i in range(D+1): s[i]+=g[F][i]
    return s,len(ne)
H=[1,1,3,4,4,5,6,7]
prod=[0]*(D+1); prod[0]=1
for h in H: prod=mul(prod,geom(h))
for nm,less in (('A  earlier-lower',A),('B  later-lower (DKK)',B)):
    s,nf=rpp_series(less)
    print(f"{nm}: {nf} nonempty filters; series {s[:12]} ...  matches hook product: {s==prod}")
print(f"hook product for {H}: {prod[:12]} ...")
