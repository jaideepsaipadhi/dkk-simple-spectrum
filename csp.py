"""
GLOBAL construction: stop sweeping the word, solve for the poset directly.

Motivation: in Ex 5.12 the element 1@1 sits below 2@2 but is created by a later
word step, so no sequential sweep can produce it.  Instead search for level
assignments satisfying the structural constraints, using v (known exactly from
the c-trace) as the only input besides the diagram and framing vertex.

Constraints:
  C1  |L_i| = v_i
  C2  distinct levels within a colour
  C3  level 0 holds exactly one element, of the framing colour k
  C4  every element at level d>0 has a neighbour-coloured element at level d-1
  C5  (Dranowski Prop 2.7) between CONSECUTIVE elements of a colour fibre lie
      exactly two neighbours, in the poset order induced by the edge set
      (edges = neighbouring colour, consecutive level)
"""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from trueheap import true_heap, is_lambda_min


def solutions(rs, v, k, maxlev=14, cap=20000):
    """DFS over rows: row d is a set of colours."""
    out = []
    rows = []
    remaining = dict(v)
    seen = [0]

    def rec(d, prev):
        seen[0] += 1
        if seen[0] > cap:
            return
        if all(x == 0 for x in remaining.values()):
            out.append([set(r) for r in rows])
            return
        if d > maxlev:
            return
        avail = [i for i in rs.I if remaining[i] > 0 and
                 (d == 0 or any(j in prev for j in rs.adj[i]))]
        # choose a nonempty subset of avail for this row
        n = len(avail)
        for mask in range(1, 1 << n):
            S = {avail[t] for t in range(n) if mask >> t & 1}
            if d == 0 and S != {k}:
                continue
            for i in S:
                remaining[i] -= 1
            rows.append(S)
            rec(d + 1, S)
            rows.pop()
            for i in S:
                remaining[i] += 1
        if out:
            return

    rec(0, set())
    return out


def rows_to_levels(rows):
    L = {}
    for d, S in enumerate(rows):
        for i in S:
            L.setdefault(i, []).append(d)
    return {i: sorted(x) for i, x in L.items()}


def prop27_ok(rs, L):
    """Between consecutive fibre elements, exactly two neighbours strictly
    between, counting only levels in the open interval."""
    for i, ls in L.items():
        for a, b in zip(ls, ls[1:]):
            mids = sum(1 for j in rs.adj[i] for l in L.get(j, [])
                       if a < l < b)
            if mids != 2:
                return False
    return True


if __name__ == "__main__":
    tests = [
        ("Ex 5.10", ("D",4,2), (1,2,3,4,2)),
        ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2)),
        ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2)),
    ]
    for label, (typ,n,k), word in tests:
        rs = RootSystem(typ, n); lam = rs.fundamental(k)
        v = dict(zip(rs.I, v_from_trace(rs, trace(rs, word, lam))))
        sols = solutions(rs, v, k)
        good = [rows_to_levels(s) for s in sols]
        good = [g for g in good if prop27_ok(rs, g)]
        print(f"{label}: v={tuple(v[i] for i in rs.I)}  "
              f"{len(sols)} level-assignments satisfy C1-C4, "
              f"{len(good)} also satisfy Prop 2.7")
        for g in good[:4]:
            print(f"     {dict(sorted(g.items()))}")
