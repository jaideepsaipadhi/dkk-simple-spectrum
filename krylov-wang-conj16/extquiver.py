"""
Extended quivers of trivial type (Krylov-Wang, Section 3.4), generated exactly
from the recursive definition, so that we can compute

        { Hesselink-type }  \  { V(A) : A of trivial type }

for G = SL2^k.  Section 3.4 conjectures that for a product of copies of sl2
the representations satisfying (2) are exactly the V(A) of trivial type, so the
set difference is the natural place for a counterexample to Conjecture 1.6.

An extended quiver is A = (X, Y) with X a finite vertex set and Y a multiset of
  * elements a of X            -- a framing at a,          summand V_1 at a
  * unordered pairs {a,b}      -- an edge (a=b: self-loop), V_2 at a if a=b,
                                  else V_1 box V_1 at a,b
  * unordered triples {a,b,c}  -- a 3-hyperedge,           V_1 box V_1 box V_1

Trivial type, from the paper:
  1. a single vertex, a vertex with a framing, a path graph with a single
     framing on each terminal vertex (a path of length 1 has two framings on
     the sole vertex), an l-cycle (a 1-cycle is a self-loop), or a 3-hyperedge;
  2. an extended quiver of trivial type with a leaf attached to one of its
     vertices;
  3. a disjoint union of extended quivers of trivial type.
"""
import itertools
from functools import lru_cache


def hw_frame(k, a):
    z = [0] * k; z[a] = 1; return tuple(z)


def hw_loop(k, a):
    z = [0] * k; z[a] = 2; return tuple(z)


def hw_edge(k, a, b):
    z = [0] * k; z[a] = 1; z[b] = 1; return tuple(z)


def hw_tri(k, a, b, c):
    z = [0] * k; z[a] = z[b] = z[c] = 1; return tuple(z)


def basic_connected(k, S):
    """Rule 1: basic trivial-type quivers on exactly the vertex set S."""
    S = tuple(sorted(S))
    out = []
    if len(S) == 1:
        a = S[0]
        out.append([])                                   # single vertex
        out.append([hw_frame(k, a)])                     # framed vertex
        out.append([hw_frame(k, a), hw_frame(k, a)])     # path of length 1
        out.append([hw_loop(k, a)])                      # 1-cycle (self-loop)
        return out
    # paths: any ordering of S (up to reversal), framing on both ends
    for perm in itertools.permutations(S):
        if perm[0] > perm[-1]:
            continue
        Y = [hw_edge(k, perm[i], perm[i + 1]) for i in range(len(perm) - 1)]
        Y += [hw_frame(k, perm[0]), hw_frame(k, perm[-1])]
        out.append(Y)
    # l-cycles: cyclic orderings of S
    if len(S) == 2:
        a, b = S
        out.append([hw_edge(k, a, b), hw_edge(k, a, b)])          # 2-cycle
    else:
        first = S[0]
        for perm in itertools.permutations(S[1:]):
            cyc = (first,) + perm
            if len(cyc) > 2 and cyc[1] > cyc[-1]:
                continue
            Y = [hw_edge(k, cyc[i], cyc[(i + 1) % len(cyc)])
                 for i in range(len(cyc))]
            out.append(Y)
    # 3-hyperedge
    if len(S) == 3:
        out.append([hw_tri(k, *S)])
    return out


def connected_trivial(k, S):
    """Rules 1 and 2: connected trivial-type quivers on the vertex set S."""
    S = tuple(sorted(S))
    out = list(basic_connected(k, S))
    # rule 2: a trivial-type quiver on a proper subset, plus leaves
    for r in range(1, len(S)):
        for sub in itertools.combinations(S, r):
            rest = [v for v in S if v not in sub]
            for base in connected_trivial(k, sub):
                # attach the remaining vertices as leaves, in any order,
                # each to any already-present vertex
                def attach(cur, present, remaining):
                    if not remaining:
                        out.append(list(cur))
                        return
                    for i, v in enumerate(remaining):
                        for anchor in present:
                            attach(cur + [hw_edge(k, anchor, v)],
                                   present + [v],
                                   remaining[:i] + remaining[i + 1:])
                attach(base, list(sub), rest)
    return out


def trivial_type_reps(k):
    """Rule 3 as well: all V(A) for A of trivial type on k labelled vertices,
    as sorted tuples of highest weights."""
    seen = set()
    verts = list(range(k))

    def partitions(items):
        if not items:
            yield []
            return
        first, rest = items[0], items[1:]
        for r in range(len(rest) + 1):
            for chosen in itertools.combinations(rest, r):
                block = [first] + list(chosen)
                remain = [x for x in rest if x not in chosen]
                for tail in partitions(remain):
                    yield [block] + tail

    for part in partitions(verts):
        pieces = [connected_trivial(k, tuple(b)) for b in part]
        for combo in itertools.product(*pieces):
            Y = []
            for c in combo:
                Y += c
            if Y:
                seen.add(tuple(sorted(Y)))
    return seen


if __name__ == "__main__":
    for k in (1, 2, 3):
        tt = trivial_type_reps(k)
        print(f"k = {k}: {len(tt)} distinct V(A) of trivial type")
        if k <= 2:
            for x in sorted(tt):
                print("   ", list(x))
