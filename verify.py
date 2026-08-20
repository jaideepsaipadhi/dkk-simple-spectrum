"""
verify.py -- one-command verification of every headline claim.

    python3 verify.py

Checks, in order:
  1. the paper's Examples 5.10-5.12 (mu, v, and the printed posets)
  2. the socle construction reproduces those posets
  3. the criterion agrees with the independent quiver-variety solver
  4. well-definedness (independent of the reduced word)
  5. the type-D theorem  simple <=> r + 2m <= 2
"""
import time, sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from proj2 import Proj2
from gls import iterated_socle
from dfs3 import solve
from weylprops import reduced_words
from epsilon import eps_of_mu

def socv(rs, k):
    pr = Proj2(rs, k)
    return max((i for (d, i), m in pr.dim.items() if d == pr.top and m), default=None)

def setup(typ, n, k):
    rs = RootSystem(typ, n)
    nu = {j: socv(rs, j) for j in rs.I}
    inv = {b: a for a, b in nu.items()}
    return rs, Proj2(rs, inv[k]), [i for i in rs.I if len(rs.adj[i]) >= 3][0]

def poset(rs, pr, word, lam):
    tr = trace(rs, word, lam)
    seq = []
    for _, i, c, _ in tr:
        seq += [i] * c
    X = iterated_socle(pr, rs, seq)
    rows = {}
    for d in range(pr.top + 1):
        for i in rs.I:
            m = len(X.get((d, i), []))
            if m:
                rows.setdefault(d, []).extend([i] * m)
    return list(reversed([sorted(rows[d]) for d in sorted(rows)])), X

def criterion(rs, pr, t, word, lam):
    _, X = poset(rs, pr, word, lam)
    return all(len(X.get((d, t), [])) <= 1 for d in range(pr.top + 1))

PAPER = [("Ex 5.10", ("D",4,2), (1,2,3,4,2), (2,2,1,1), [[2],[1,3,4],[2],[1]]),
         ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2), (1,3,2,1,1),
          [[2],[1,3],[2,4,5],[3],[2]]),
         ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2), (2,3,2,1,1),
          [[2],[1,3],[2,4,5],[3],[2],[1]])]

ok = True
print("1. DKK Examples 5.10-5.12: dimension vectors")
for lab, (typ,n,k), word, vexp, rows in PAPER:
    rs = RootSystem(typ,n)
    v = tuple(int(x) for x in v_from_trace(rs, trace(rs, word, rs.fundamental(k))))
    good = v == vexp
    ok &= good
    print(f"   {lab}: v = {v}  expected {vexp}   {'OK' if good else 'FAIL'}")

print("\n2. The socle construction reproduces the printed posets")
for lab, (typ,n,k), word, vexp, rows in PAPER:
    rs, pr, t = setup(typ,n,k)
    got, _ = poset(rs, pr, word, rs.fundamental(k))
    good = got == rows
    ok &= good
    print(f"   {lab}: {got}   {'OK' if good else 'FAIL'}")

print("\n3. Criterion vs the independent solver")
t0 = time.time()
for typ, n, k in [("D",4,2), ("D",5,2), ("D",6,3), ("E",6,3), ("E",6,2)]:
    rs, pr, t = setup(typ,n,k); lam = rs.fundamental(k)
    agree = tot = 0
    for mu, word in orbit_with_words(rs, lam).items():
        if not word: continue
        v = tuple(int(x) for x in v_from_trace(rs, trace(rs, word, lam)))
        r = solve(rs, k, v)
        if not r: continue
        g = list(r.values())[0]
        truth = all(c <= 1 for row in g.values() for c in row.values())
        tot += 1
        if criterion(rs, pr, t, word, lam) == truth: agree += 1
    ok &= (agree == tot)
    print(f"   {rs.name()} om{k}: {agree}/{tot}   {'OK' if agree==tot else 'FAIL'}")
print(f"   ({time.time()-t0:.0f}s)")

print("\n4. Well-definedness (independent of the reduced word)")
for typ, n, k in [("D",4,2), ("D",5,2)]:
    rs, pr, t = setup(typ,n,k); lam = rs.fundamental(k)
    agree = tot = 0
    for mu, word in orbit_with_words(rs, lam).items():
        if not word or len(word) > 7: continue
        rw = reduced_words(rs, word)
        if len(rw) > 30: continue
        vals = {criterion(rs, pr, t, w2, lam) for w2 in rw}
        tot += 1; agree += (len(vals) == 1)
    ok &= (agree == tot)
    print(f"   {rs.name()} om{k}: {agree}/{tot} elements word-independent   "
          f"{'OK' if agree==tot else 'FAIL'}")

print("\n5. Type-D theorem:  simple  <=>  r + 2m <= 2")
for n in (4,5,6):
    rs = RootSystem("D",n); t = n-2
    agree = tot = 0
    for k in range(1, n-1):
        lam = rs.fundamental(k)
        for mu, word in orbit_with_words(rs, lam).items():
            if not word: continue
            v = [int(x) for x in v_from_trace(rs, trace(rs, word, lam))]
            e = eps_of_mu(n, k, v)
            r_ = sum(1 for i in (n-2, n-1) if e[i] != 0)
            m_ = sum(1 for i in range(0, n-2) if e[i] == -1)
            res = solve(rs, k, tuple(v))
            if not res: continue
            g = list(res.values())[0]
            truth = all(c <= 1 for row in g.values() for c in row.values())
            tot += 1
            agree += (truth == (r_ + 2*m_ <= 2))
    ok &= (agree == tot)
    print(f"   D{n}: {agree}/{tot}   {'OK' if agree==tot else 'FAIL'}")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
sys.exit(0 if ok else 1)
