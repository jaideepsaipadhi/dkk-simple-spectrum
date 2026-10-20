"""Is the counterexample robust to the orientation convention?

The true hook multiset is computed where the RPP series IS a hook product
(orientation A, earlier-letter-lower).  Def 2.6's recursion is stated in the
corner-as-maximum convention, so the project applies it to the dual.  If the
true multiset is unreachable in BOTH orientations, the conclusion does not
depend on that choice.
usage: python3 orient.py
"""
import sys, itertools, collections
sys.path.insert(0,'.')
from nsl import cartan, heap_down
from novel import dual, relabel_natural
from no26 import leq_matrix, dk_intervals
from filterhook import hooks_of
from choice import hooks_from_choice
W=[1,0,1,2,3,1,0,2]
down=heap_down(cartan('D',4),W)
A=relabel_natural(down)             # as built: earlier letter lower
Bo=relabel_natural(dual(down))      # the dual
true=hooks_of(A)
print(f"true hook multiset (orientation A, verified independently): {true}")
for nm,P in (('A  as built',A),('B  dual',Bo)):
    leq=leq_matrix(P)
    tops=collections.defaultdict(list)
    for (v,u,x,y) in dk_intervals(P,leq): tops[u].append((v,x,y))
    keys=sorted(tops); nc=1
    for u in keys: nc*=(len(tops[u])+1)
    got=set()
    for pick in itertools.product(*[[None]+list(range(len(tops[u]))) for u in keys]):
        ch={u:tops[u][p] for u,p in zip(keys,pick) if p is not None}
        got.add(tuple(sorted(hooks_from_choice(P,leq,ch))))
    hit=tuple(true) in got
    print(f"\n{nm}: elements topping a d_k-interval {keys}, {nc} choices")
    for g in sorted(got): print(f"    reachable multiset {list(g)}")
    print(f"    true multiset reachable: {hit}")
