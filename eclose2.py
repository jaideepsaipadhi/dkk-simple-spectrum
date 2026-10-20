"""The 14 spread cases against the non-simple weights sharing their coarse
invariants.  v determines mu determines w, so a criterion must be a function
of v -- the point is to find a SIMPLE one.
usage: python3 eclose2.py
"""
import collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
rs=RootSystem('E',6); t=4
cand=[]
for k in (3,4,5):
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        if dv[t-1]!=3: continue
        simple=all(c<=1 for r in rows.values() for c in r.values())
        tr=trace(rs,w,lam)
        # degrees in sigma at which the node steps occur (position in the word)
        pos=[j for j,(_,i,c,_) in enumerate(tr) if i==t]
        cand.append(dict(ka=k,simple=simple,v=dv,mu=tuple(mu),L=len(w),
                         pos=tuple(pos),word=w,
                         seq=tuple(i for _,i,c,_ in tr)))
S=[r for r in cand if r['simple']]; NS=[r for r in cand if not r['simple']]
print(f"{len(S)} spread/simple, {len(NS)} fill/not-simple\n")
print("the 14, by kappa:")
for r in S: print(f"  k={r['ka']} L={r['L']:2} v={r['v']} nodepos={r['pos']}")
print("\n-- candidate coarse invariants --")
for key,name in [(lambda r:r['pos'],'positions of the node steps in sigma'),
                 (lambda r:tuple(b-a for a,b in zip(r['pos'],r['pos'][1:])),'gaps between node steps'),
                 (lambda r:r['L'],'word length'),
                 (lambda r:(r['ka'],r['L']),'(kappa, length)'),
                 (lambda r:max(r['v']),'max_i v_i'),
                 (lambda r:sum(r['v']),'sum v_i')]:
    d=collections.defaultdict(lambda:[0,0])
    for r in cand: d[key(r)][r['simple']]+=1
    mixed=[kk for kk,vv in d.items() if vv[0] and vv[1]]
    print(f"\n  {name}: {len(d)} values, {len(mixed)} mixed"
          f" -> {'*** SEPARATES ***' if not mixed else 'mixed'}")
    if not mixed or len(d)<=12:
        for kk in sorted(d,key=str): print(f"     {kk!s:28} not={d[kk][0]:3} simple={d[kk][1]:3}")
