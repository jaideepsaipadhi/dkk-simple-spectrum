"""Does SOME choice of d_k-interval rescue the skew formula on a non-d-complete witness?

Outside the d-complete class an element may top several d_k-intervals.  Def 2.6's hook
recursion and Def 3.1's excitation both then depend on which one is used -- a choice the
earlier code made arbitrarily (tops[u][0]).  Here we enumerate EVERY choice function,
recompute both the hook function and the excited-diagram graph from it, and test the skew
identity on all order filters.
"""
import sys, itertools, collections
sys.path.insert(0, '.')
from nsl import cartan, heap_down
from novel import relabel_natural, covers
from no26 import leq_matrix, dk_intervals
from filterhook import hooks_of
from skew import to_proctor, series_prod
from skew2 import heap_sub_series

N = 18

def hooks_from_choice(dP, leq, choice):
    n = len(dP)
    order = sorted(range(n), key=lambda u: sum(1 for w in range(n) if leq[w][u]))
    h = {}
    for u in order:
        if u in choice:
            v, x, y = choice[u]
            h[u] = h[x] + h[y] - h[v]
        else:
            h[u] = sum(1 for w in range(n) if leq[w][u])
    return [h[u] for u in range(n)]

def excited_from_choice(dP, leq, choice, F, cov, cap=30000):
    n = len(dP)
    moves = {}
    for u, (v, x, y) in choice.items():
        I = [w for w in range(n) if leq[v][w] and leq[w][u]]
        rim = set(w for w in I if w in cov[u]) | set(w for w in I if u != w and v in cov[w])
        moves[u] = (v, frozenset(rim))
    start = (frozenset(F), frozenset()); seen = {start}; st = [start]
    while st:
        D, B = st.pop()
        for u in list(D):
            if u not in moves: continue
            v, rim = moves[u]
            if v in D or any(z in D for z in rim): continue
            s = (frozenset((set(D) - {u}) | {v}), frozenset((set(B) - set(rim)) | {u}))
            if s not in seen:
                seen.add(s); st.append(s)
                if len(seen) > cap: return None
    return seen

def run(typ, rk, w):
    C = cartan(typ, rk)
    down = heap_down(C, list(w)); n = len(down)
    dn = relabel_natural(down); th = hooks_of(dn)
    dP, col, order = to_proctor(down, list(w))
    leq = leq_matrix(dP)
    cov = [[v for v in range(n) if v != u and leq[v][u] and
            not any(z not in (u, v) and leq[v][z] and leq[z][u] for z in range(n))]
           for u in range(n)]
    tops = collections.defaultdict(list)
    for (v, u, x, y) in dk_intervals(dP, leq): tops[u].append((v, x, y))
    amb = {u: v for u, v in tops.items()}
    nchoice = 1
    for u in amb: nchoice *= (len(amb[u]) + 1)      # +1 = "use no interval at u"
    print(f"{typ}{rk} w={''.join('s'+str(i) for i in w)} |P|={n} true hooks={th}")
    print(f"   elements topping a d_k-interval: "
          f"{ {u: len(v) for u, v in amb.items()} };  {nchoice} choice functions")
    ups = [frozenset(x for x in range(n) if leq[u][x]) for u in range(n)]
    filters = sorted({frozenset().union(*[ups[u] for u in S])
                      for r in (1, 2) for S in itertools.combinations(range(n), r)},
                     key=lambda f: len(f))
    filters = [f for f in filters if 0 < len(f) < n]
    lhs = {F: heap_sub_series(down, sorted(order[u] for u in range(n) if u not in F), N)
           for F in filters}
    keys = sorted(amb)
    good = 0; best = -1
    for pick in itertools.product(*[[None] + list(range(len(amb[u]))) for u in keys]):
        choice = {u: amb[u][p] for u, p in zip(keys, pick) if p is not None}
        h = hooks_from_choice(dP, leq, choice)
        if sorted(h) != th: continue                  # F = empty forces the true multiset
        okc = 0
        for F in filters:
            E = excited_from_choice(dP, leq, choice, F, cov)
            if E is None: break
            rhs = [0]*(N+1)
            for D, B in E:
                sh = sum(h[v] for v in B)
                if sh > N: continue
                pr = series_prod([h[u] for u in range(n) if u not in D], N)
                for d in range(N+1-sh): rhs[d+sh] += pr[d]
            if rhs[:N+1] == lhs[F][:N+1]: okc += 1
        best = max(best, okc)
        if okc == len(filters): good += 1
    print(f"   choices giving the TRUE hook multiset and satisfying ALL {len(filters)} "
          f"filters: {good}")
    print(f"   best any choice achieves: {best} of {len(filters)} filters\n")

if __name__ == "__main__":
    run('D', 4, (1, 0, 1, 2, 3, 1, 0, 2))
    run('D', 4, (0, 1, 2, 1, 0, 3, 1, 0, 2))
