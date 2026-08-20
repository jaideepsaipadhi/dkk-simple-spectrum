"""
Constraint harvest.

Enumerate the Weyl group, compute the reduced-word c-trace against a
fundamental weight lambda, and tabulate:
  - the c-multiset (how far w is from lambda-minuscule)
  - whether the c-multiset depends on the choice of reduced word
  - correlation with fully commutative / minuscule
"""

from collections import deque, Counter
from rootsys import RootSystem
from wordtrace import trace, v_from_trace


def enumerate_W(rs, cap=None):
    """BFS over W. Element key = images of all fundamental weights.
    Returns dict key -> (reduced_word, length)."""
    idw = tuple(tuple(rs.fundamental(i)) for i in rs.I)
    seen = {idw: ()}
    q = deque([idw])
    while q:
        cur = q.popleft()
        if cap is not None and len(seen[cur]) >= cap:
            continue
        for i in rs.I:
            nxt = tuple(tuple(rs.act(i, mu)) for mu in cur)
            if nxt not in seen:
                seen[nxt] = (i,) + seen[cur]
                q.append(nxt)
    return seen


def all_reduced_words(rs, word):
    """All reduced words, via braid + commutation moves."""
    start = tuple(word)
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for k in range(len(u) - 1):
            a, b = u[k], u[k + 1]
            if a == b:
                continue
            if b not in rs.adj[a]:
                nxt = u[:k] + (b, a) + u[k + 2:]
                if nxt not in seen:
                    seen.add(nxt); q.append(nxt)
            elif k + 2 <= len(u) - 1 and u[k + 2] == a:
                nxt = u[:k] + (b, a, b) + u[k + 3:]
                if nxt not in seen:
                    seen.add(nxt); q.append(nxt)
    return seen


def cprofile(rs, word, lam):
    tr = trace(rs, word, lam)
    return tuple(sorted((c for _, _, c, _ in tr), reverse=True)), v_from_trace(rs, tr)


def study(typ, n, k, wordcap=None, check_indep=True, sample=400):
    rs = RootSystem(typ, n)
    lam = rs.fundamental(k)
    W = enumerate_W(rs, cap=wordcap)
    print(f"=== {rs.name()}  lambda = omega_{k}   |W enumerated| = {len(W)}")

    prof_count = Counter()
    maxc_count = Counter()
    examples = {}
    indep_fail = []

    items = list(W.items())
    for idx, (key, word) in enumerate(items):
        prof, v = cprofile(rs, word, lam)
        maxc = prof[0] if prof else 1
        prof_count[prof] += 1
        maxc_count[maxc] += 1
        if maxc not in examples:
            examples[maxc] = (word, v, prof)
        # independence of reduced word (sampled -- the move graph can be large)
        if check_indep and idx < sample and 0 < len(word) <= 9:
            rws = all_reduced_words(rs, word)
            if len(rws) <= 200:
                profs = {cprofile(rs, u, lam)[0] for u in rws}
                if len(profs) > 1:
                    indep_fail.append((word, profs))

    print(f"  max c distribution: "
          + ", ".join(f"c_max={m}: {cnt}" for m, cnt in sorted(maxc_count.items())))
    for m in sorted(examples):
        wd, v, prof = examples[m]
        s = ''.join('s' + str(i) for i in wd)
        print(f"    c_max={m}  e.g. w={s or 'e'}  v={v}  c-profile={prof}")
    if check_indep:
        if indep_fail:
            print(f"  !! c-profile DEPENDS on reduced word for "
                  f"{len(indep_fail)} sampled elements, e.g.:")
            for wd, profs in indep_fail[:3]:
                print(f"       w={''.join('s'+str(i) for i in wd)}  profiles={profs}")
        else:
            print("  c-profile independent of reduced word on all sampled elements")
    print()
    return maxc_count


if __name__ == "__main__":
    print("Constraint harvest: c-profiles across Weyl groups\n")
    study("D", 4, 2)
    study("D", 5, 2)
    study("D", 4, 1)
    study("A", 4, 2)
