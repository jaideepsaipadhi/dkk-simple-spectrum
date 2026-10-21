"""End-to-end validation of Table 1 of appendixE: every row recomputed from the
printed word and compared with the printed data.
usage: python3 validate_nine.py
"""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
from weylprops import is_fully_commutative, reduced_words
ROWS=[
 (1,'E6',4,[2,4,5,6,3,4,5,1,3,4],10,[1]*9+[3],(0,-3,0,2,0,0),(1,3,2,3,2,1),12),
 (2,'E7',3,[3,4,5,6,7,2,4,5,6,1,3,4,5,2,4,3],16,[1]*6+[2]+[1]*8+[3],(1,0,-3,2,0,0,0),(2,2,5,4,3,2,1),19),
 (3,'E7',3,[1,3,4,5,6,7,2,4,5,6,3,4,5,2,4,3],16,[1]*15+[3],(-3,0,2,0,0,0,0),(3,2,3,4,3,2,1),18),
 (4,'E7',3,[1,3,4,5,6,7,2,4,5,6,1,3,4,5,2,4,3],17,[1]*6+[2]+[1]*8+[3,1],(-1,0,-2,2,0,0,0),(3,2,5,4,3,2,1),20),
 (5,'E7',4,[2,4,5,6,3,4,5,1,3,4],10,[1]*9+[3],(0,-3,0,2,0,0,1),(1,3,2,3,2,1,0),12),
 (6,'E7',4,[7,2,4,5,6,3,4,5,1,3,4],11,[1]*9+[3,1],(0,-3,0,2,0,1,-1),(1,3,2,3,2,1,1),13),
 (7,'E7',4,[6,7,2,4,5,6,3,4,5,1,3,4],12,[1]*9+[3,1,1],(0,-3,0,2,1,-1,0),(1,3,2,3,2,2,1),14),
 (8,'E7',4,[5,6,7,2,4,5,6,3,4,5,1,3,4],13,[1]*9+[3,1,1,1],(0,-3,0,3,-1,0,0),(1,3,2,3,3,2,1),15),
 (9,'E7',5,[2,4,5,6,7,3,4,5,6,1,3,4,5],13,[1]*12+[3],(0,-3,0,2,0,0,0),(1,3,2,3,3,2,1),15),
]
bad=0
for idx,typ,k,w,L,ctr,mu_p,v_p,dim_p in ROWS:
    rs=RootSystem(typ[0],int(typ[1:])); lam=rs.fundamental(k)
    errs=[]
    if len(w)!=L: errs.append(f"length {len(w)} != {L}")
    # minimal-length rep check: w must appear in orbit_with_words for this lam
    orb=orbit_with_words(rs,lam)
    mu=tuple(rs.act_word(w,lam))
    if mu not in orb: errs.append("mu not in the orbit map")
    elif len(orb[mu])!=L: errs.append(f"not minimal: orbit gives length {len(orb[mu])}")
    if mu!=tuple(mu_p): errs.append(f"mu {mu} != printed {tuple(mu_p)}")
    tr=trace(rs,w,lam); cs=[c for _,_,c,_ in tr]
    if cs!=ctr: errs.append(f"c-trace {cs} != printed {ctr}")
    _,rows,dv=soc_sigma(rs,k,w)
    if dv!=tuple(v_p): errs.append(f"v {dv} != printed {tuple(v_p)}")
    if sum(dv)!=dim_p: errs.append(f"dim {sum(dv)} != printed {dim_p}")
    simple=all(c<=1 for r in rows.values() for c in r.values())
    if not simple: errs.append("NOT simple spectrum")
    if max(cs)!=3: errs.append(f"max c = {max(cs)} != 3")
    fc=is_fully_commutative(rs,w) if L<=13 else None
    print(f"({idx}) {typ} om{k}: " + ("OK" if not errs else "FAIL " + "; ".join(errs))
          + (f"   fully commutative: {fc}" if fc is not None else "   (fc not tested, long)"))
    bad+=len(errs)
print(f"\nrow errors: {bad}")
