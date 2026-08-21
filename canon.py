"""
THE REPRESENTATIVES ARE CANONICAL: a construction of V(mu_a) with no word.

reduce.py cuts (*o) down to the k weights mu_a whose only descent is t-1, and
observes that at those weights V(mu_a) is the greedy top truncation of I_kappa
at every vertex -- statement (B).  This script identifies the module that the
truncation must be, intrinsically.

A graded top truncation is not determined by its multiplicities alone: where
I_kappa carries multiplicity 2 and the truncation takes only 1, a LINE has to be
chosen.  Over all extremal modules exactly three lines ever occur in such a
slot, and the representatives always use the same one.  It is the canonical one:
tau, the diagram automorphism exchanging the two leaves, fixes I_kappa (kappa is
a chain vertex), acts on each multiplicity-2 slot as an involution with
eigenvalues +1 and -1, and the line the representatives use is the +1
eigenline.  The other two lines are the two "special" ones, kernels of the two
leaf arrows at t, which tau exchanges.

So define, for a dimension vector v, the CANONICAL TOP TRUNCATION T(v): at each
vertex j take I_kappa[d] in full for every degree above the threshold, nothing
below it, and at the threshold degree -- if only part of a multiplicity-2 slot
is wanted -- the tau-fixed line.  Then

    (B)   is equivalent to    T(v_a) is a submodule of I_kappa,

because dim T(v_a) = v_a by construction, and a submodule with a w-extremal
dimension vector is unique (Baumann-Kamnitzer-Tingley), hence equal to V(mu_a).
This replaces "run the socle construction along a reduced word and look at the
answer" by "check that an explicitly written subspace is closed under the
arrows".

tau is realised on Pi e_{kappa'} as follows.  Proj2 builds P_{d+1}[j] as a
quotient of (+)_{i ~ j} P_d[i], and the sign convention eps(a,b) = +-1 is
tau-symmetric because n-1 and n are both greater than t.  So the construction is
tau-equivariant, and tau is determined on P_{d+1}[j] by
tau(a_{i->j} x) = a_{tau i -> tau j}(tau x), the arrow images spanning.

Closure of T(v) is then almost formal, because of two facts about I_kappa and
tau alone -- no words, no dimension vectors -- both verified here without
exception:

  (T1)  if the image of a full slot in a multiplicity-2 chain slot is only
        1-dimensional, that image IS the tau-line                  (22/22);
  (T2)  the tau-line of a multiplicity-2 chain slot is annihilated by exactly
        the arrows into 1-dimensional CHAIN slots -- it is nonzero into
        2-dimensional chain slots (36/36) and into both leaves (44/44), and
        zero into 1-dimensional chain slots (22/22).

Every closure obligation for T(v_a) is then of one of four kinds:

    source full, target full          nothing to check          (437)
    source full, target partial       the image is the tau-line, by (T1)  (10)
    source the tau-line, target full  nothing to check           (72)
    source the tau-line, target partial   the image is tau-fixed, hence in
        the target tau-line, by equivariance                      (18)
    source the tau-line, target empty     the image vanishes, by (T2)     (12)

so what is left of (B) is the arithmetic that the greedy thresholds of v_a never
produce any other pattern -- in particular never send a full slot into an empty
one.  This script verifies that too, by classifying every obligation.

This script builds tau, checks it is an involutive automorphism, verifies (T1)
and (T2), and verifies that T(v_a) is a submodule equal to V(mu_a) at every
representative in range, with every closure obligation of a resolved kind.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from homological import setup
from gls import socle_step
from wordtrace import trace
from qvar import rref, P

CASES = [("D", 5, 2), ("D", 6, 2), ("D", 7, 2), ("D", 6, 3), ("D", 7, 3),
         ("D", 8, 3), ("D", 7, 4), ("D", 8, 4), ("D", 8, 5), ("D", 9, 4)]


def solve(rows, target, m):
    """Express target in the span of rows; return coefficients or None."""
    aug = [rows[i][:] + [1 if x == i else 0 for x in range(len(rows))]
           for i in range(len(rows))]
    R, piv = rref(aug, m + len(rows))
    x = target[:] + [0] * len(rows)
    for r, c in zip(R, piv):
        if c < m and x[c] % P:
            g = x[c]
            x = [(x[i] - g * r[i]) % P for i in range(m + len(rows))]
    if any(y % P for y in x[:m]):
        return None
    return [(-y) % P for y in x[m:]]


def build_tau(pr, rs, n, kk):
    """T[(d,i)] : P_d[i] -> P_d[tau i], as a list of rows."""
    sw = {n - 1: n, n: n - 1}
    tv = lambda i: sw.get(i, i)
    T = {(0, kk): [[1]]}
    for d in range(pr.top):
        for j in rs.I:
            m = pr.dim.get((d + 1, j), 0)
            if not m:
                continue
            span, img = [], []
            for i in rs.adj[j]:
                if (d, i, j) not in pr.arrow or (d, i) not in T:
                    continue
                mi = pr.dim.get((d, i), 0)
                for b in range(mi):
                    e = [1 if x == b else 0 for x in range(mi)]
                    span.append(pr.act(d, i, e, j))
                    te = T[(d, i)][b] if (d, i) in T else None
                    if te is None:
                        return None
                    img.append(pr.act(d, tv(i), te, tv(j)))
            M = []
            for b in range(m):
                tgt = [1 if x == b else 0 for x in range(m)]
                co = solve(span, tgt, m)
                if co is None:
                    return None
                row = [0] * pr.dim.get((d + 1, tv(j)), 0)
                for c, vv in zip(co, img):
                    if c % P:
                        for x in range(len(row)):
                            row[x] = (row[x] + c * vv[x]) % P
                M.append(row)
            T[(d + 1, j)] = M
    return T


def apply_tau(T, d, j, vec):
    M = T[(d, j)]
    w = [0] * (len(M[0]) if M else 0)
    for b, a in enumerate(vec):
        if a % P:
            for x in range(len(w)):
                w[x] = (w[x] + a * M[b][x]) % P
    return w


def tau_fixed(pr, rs, T, d, j, n):
    """Basis of the tau-fixed subspace of P_d[j] (j a tau-fixed vertex)."""
    m = pr.dim[(d, j)]
    eqs = []
    for b in range(m):
        e = [1 if x == b else 0 for x in range(m)]
        w = apply_tau(T, d, j, e)
        eqs.append([(w[x] - e[x]) % P for x in range(m)])
    rows = [[eqs[b][x] for b in range(m)] for x in range(m)]
    R, piv = rref(rows, m)
    free = [c for c in range(m) if c not in piv]
    out = []
    for f in free:
        x = [0] * m
        x[f] = 1
        for r, c in zip(R, piv):
            x[c] = (-r[f]) % P
        out.append(x)
    return out


def canonical(pr, rs, T, v, n):
    """T(v): greedy top truncation, tau-fixed line at a partial 2-slot."""
    N = {}
    for j in rs.I:
        rem = v[j - 1]
        for d in range(pr.top, -1, -1):
            m = pr.dim.get((d, j), 0)
            if not m or rem <= 0:
                continue
            take = min(m, rem)
            rem -= take
            if take == m:
                N[(d, j)] = [[1 if x == b else 0 for x in range(m)]
                             for b in range(m)]
            else:
                F = tau_fixed(pr, rs, T, d, j, n)
                if len(F) != take:
                    return None
                N[(d, j)] = rref(F, m)[0]
        if rem:
            return None
    return N


def is_submodule(pr, rs, N):
    for (d, i), vecs in N.items():
        for j in rs.adj[i]:
            if (d, i, j) not in pr.arrow:
                continue
            m = pr.dim.get((d + 1, j), 0)
            if not m:
                continue
            R, piv = rref([r[:] for r in N.get((d + 1, j), [])], m)
            for vec in vecs:
                x = pr.act(d, i, vec, j)
                for r, c in zip(R, piv):
                    if x[c] % P:
                        g = x[c]
                        x = [(x[y] - g * r[y]) % P for y in range(m)]
                if any(y % P for y in x):
                    return False
    return True


def module(rs, pr, word, lam):
    X = {}
    for _, i, c, _ in trace(rs, word, lam):
        for _ in range(c):
            X = socle_step(pr, rs, X, i)
    return X


def same(pr, rs, A, B):
    for d in range(pr.top + 1):
        for j in rs.I:
            m = pr.dim.get((d, j), 0)
            if not m:
                continue
            RA = rref([r[:] for r in A.get((d, j), [])], m)[0]
            RB = rref([r[:] for r in B.get((d, j), [])], m)[0]
            if RA != RB:
                return False
    return True


def norm(r):
    piv = next(i for i, y in enumerate(r) if y % P)
    inv = pow(r[piv], P - 2, P)
    return tuple((x * inv) % P for x in r)


def local_facts(pr, rs, T, n, st):
    """(T1) and (T2)."""
    for (d, i), mi in sorted(pr.dim.items()):
        if i in (n - 1, n):
            continue
        for j in rs.adj[i]:
            if (d, i, j) not in pr.arrow or not pr.dim.get((d + 1, j), 0):
                continue
            mj = pr.dim[(d + 1, j)]
            if mj == 2 and j not in (n - 1, n):
                IM = [pr.act(d, i, [1 if x == b else 0 for x in range(mi)], j)
                      for b in range(mi)]
                R = rref([r[:] for r in IM], 2)[0]
                if len(R) == 1:
                    F = tau_fixed(pr, rs, T, d + 1, j, n)
                    st['T1_tot'] += 1
                    st['T1_ok'] += (norm(R[0]) == norm(F[0]))
        if mi != 2:
            continue
        L = tau_fixed(pr, rs, T, d, i, n)
        if len(L) != 1:
            st['T2_bad'] += 1
            continue
        for j in rs.adj[i]:
            if (d, i, j) not in pr.arrow or not pr.dim.get((d + 1, j), 0):
                continue
            zero = not any(y % P for y in pr.act(d, i, L[0], j))
            want = (pr.dim[(d + 1, j)] == 1 and j not in (n - 1, n))
            st['T2_tot'] += 1
            st['T2_ok'] += (zero == want)


def obligations(pr, rs, N, st):
    """Classify every arrow obligation of the truncation N."""
    for (d, i), vecs in N.items():
        full_src = len(vecs) == pr.dim[(d, i)]
        for j in rs.adj[i]:
            if (d, i, j) not in pr.arrow:
                continue
            m = pr.dim.get((d + 1, j), 0)
            if not m:
                continue
            got = len(N.get((d + 1, j), []))
            tgt = 'full' if got == m else ('partial' if got else 'empty')
            zero = all(not any(y % P for y in pr.act(d, i, vec, j))
                       for vec in vecs)
            key = ('full' if full_src else 'line', tgt, zero)
            st['ob_' + '_'.join(map(str, key))] += 1
            st['ob_tot'] += 1
            st['ob_ok'] += key in {('full', 'full', False), ('full', 'full', True),
                                   ('full', 'partial', False),
                                   ('line', 'full', False), ('line', 'full', True),
                                   ('line', 'partial', False),
                                   ('line', 'empty', True)}


def run(typ, n, k):
    rs, pr, kk = setup(typ, n, k)
    t = n - 2
    lam = rs.fundamental(k)
    T = build_tau(pr, rs, n, kk)
    st = Counter()
    if T is None:
        st['tau_fail'] = 1
        return st
    # tau is an involutive automorphism
    for (d, j) in list(pr.dim):
        m = pr.dim[(d, j)]
        for b in range(m):
            e = [1 if x == b else 0 for x in range(m)]
            sw = {n - 1: n, n: n - 1}
            back = apply_tau(T, d, sw.get(j, j), apply_tau(T, d, j, e))
            st['inv_tot'] += 1
            st['inv_ok'] += (back == e)
            for i in rs.adj[j]:
                if (d, j, i) not in pr.arrow or not pr.dim.get((d + 1, i), 0):
                    continue
                lhs = apply_tau(T, d + 1, i, pr.act(d, j, e, i))
                rhs = pr.act(d, sw.get(j, j), apply_tau(T, d, j, e), sw.get(i, i))
                st['equi_tot'] += 1
                st['equi_ok'] += (lhs == rhs)
    local_facts(pr, rs, T, n, st)
    for mu, w in orbit_with_words(rs, lam).items():
        if mu[t - 1] < 1 or not w:
            continue
        if any(mu[j - 1] <= -1 and j not in (t - 1, t) for j in rs.I):
            continue
        X = module(rs, pr, w, lam)
        v = [sum(len(X.get((d, j), [])) for d in range(pr.top + 1)) for j in rs.I]
        N = canonical(pr, rs, T, v, n)
        st['rep'] += 1
        st['sub'] += (N is not None and is_submodule(pr, rs, N))
        st['eq'] += (N is not None and same(pr, rs, N, X))
        if N is not None:
            obligations(pr, rs, N, st)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        S += r
        print(f"  {typ}{n} k={k}: tau involutive {r['inv_ok']}/{r['inv_tot']}, "
              f"equivariant {r['equi_ok']}/{r['equi_tot']};  "
              f"T(v) a submodule {r['sub']}/{r['rep']}, "
              f"= V(mu) {r['eq']}/{r['rep']}")
    ok = (S['inv_ok'] == S['inv_tot'] and S['equi_ok'] == S['equi_tot']
          and S['sub'] == S['rep'] and S['eq'] == S['rep'] and not S['tau_fail']
          and S['T1_ok'] == S['T1_tot'] and S['T2_ok'] == S['T2_tot']
          and S['ob_ok'] == S['ob_tot'] and not S['T2_bad'])
    print(f"\n  tau involutive        {S['inv_ok']}/{S['inv_tot']}")
    print(f"  tau equivariant       {S['equi_ok']}/{S['equi_tot']}")
    print(f"  T(v_a) is a submodule {S['sub']}/{S['rep']}")
    print(f"  T(v_a) = V(mu_a)      {S['eq']}/{S['rep']}")
    print(f"  (T1)                  {S['T1_ok']}/{S['T1_tot']}")
    print(f"  (T2)                  {S['T2_ok']}/{S['T2_tot']}")
    print(f"  closure obligations   {S['ob_ok']}/{S['ob_tot']}, by kind:")
    for a, b in sorted(S.items()):
        if a.startswith('ob_') and a not in ('ob_tot', 'ob_ok'):
            print(f"      {a[3:]:24} {b}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("CANONICAL FORM VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
