"""Full commutativity without enumerating all reduced words.

w is fully commutative iff no reduced word for w contains s_i s_j s_i with
i ~ j.  If w is FC its reduced words form a single commutation class, so it
suffices to close the given word under COMMUTATION moves only and look for a
braid pattern in the class; if one appears, w is not FC, and if none does,
no braid move is ever available, so the class is everything and w is FC.
usage: python3 fctest.py
"""
from rootsys import RootSystem
def fc(rs,w):
    w=tuple(w); seen={w}; st=[w]; 
    while st:
        u=st.pop()
        for p in range(len(u)-2):
            if u[p]==u[p+2] and u[p+1] in rs.adj[u[p]]: return False,len(seen)
        for p in range(len(u)-1):
            a,b=u[p],u[p+1]
            if a!=b and b not in rs.adj[a]:
                v=u[:p]+(b,a)+u[p+2:]
                if v not in seen: seen.add(v); st.append(v)
    return True,len(seen)
CASES=[(1,'E6',[2,4,5,6,3,4,5,1,3,4]),
       (2,'E7',[3,4,5,6,7,2,4,5,6,1,3,4,5,2,4,3]),
       (3,'E7',[1,3,4,5,6,7,2,4,5,6,3,4,5,2,4,3]),
       (4,'E7',[1,3,4,5,6,7,2,4,5,6,1,3,4,5,2,4,3]),
       (5,'E7',[2,4,5,6,3,4,5,1,3,4]),
       (6,'E7',[7,2,4,5,6,3,4,5,1,3,4]),
       (7,'E7',[6,7,2,4,5,6,3,4,5,1,3,4]),
       (8,'E7',[5,6,7,2,4,5,6,3,4,5,1,3,4]),
       (9,'E7',[2,4,5,6,7,3,4,5,6,1,3,4,5])]
res={}
for idx,typ,w in CASES:
    rs=RootSystem(typ[0],int(typ[1:]))
    ok,n=fc(rs,w); res[idx]=ok
    print(f"({idx}) {typ}: fully commutative = {ok}   (commutation class size {n})")
notfc=[i for i,v in res.items() if not v]
print(f"\nnot fully commutative: {notfc}")
print(f"claim '(4) is the only one' : {notfc==[4]}")
