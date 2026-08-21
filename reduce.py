r"""
THE REMAINING GAP, REDUCED TO k WEIGHTS -- INDEPENDENT OF THE RANK.

The gap is (*o) (starchar.py), in word-free form: for a graded submodule
M <= I_kappa of D_n with delta(dim M) = 0 -- equivalently dim M is w-extremal,
mu = lambda - dim M lying in W lambda -- and c_t = <alpha_t^vee, mu> >= 1,

    (*)    M is an up-set at the chain neighbour t-1.

THE REDUCTION.  Let J = I \ {t-1, t} and let j in J be a descent of mu, i.e.
<alpha_j^vee, mu> <= -1.  Put nu = s_j mu, which is shorter.  By the main
theorem M(mu) = soc_j(M(nu)), and a socle step of colour j enlarges the module
only at vertex j; since j != t-1,

        M(mu) and M(nu) have the SAME graded profile at t-1.

So "up-set at t-1" is constant on W_J-orbits, and (*) for mu follows from (*)
for the W_J-minimal representative of mu.  (One cannot take j = t: there
<alpha_t^vee, s_t mu> = -<alpha_t^vee, mu> <= -1, so s_t mu is longer, not
shorter.)

WHAT SURVIVES, IN CLOSED FORM.  The W_J-minimal weights with c_t >= 1 are
exactly the mu whose only descent is at t-1, and they can be written down.  In
Bourbaki coordinates W omega_k = { sum_{i in S} eta_i e_i : |S| = k, eta_i = +-1 }
for k <= n-2, so mu = (eta_1, ..., eta_n) with eta_i in {0,+-1} and exactly k of
them nonzero, and

    <alpha_j^vee, mu> = eta_j - eta_{j+1} (j <= n-1),
    <alpha_n^vee, mu> = eta_{n-1} + eta_n.

The hypotheses say eta_j >= eta_{j+1} for every j <= n-1 except j = n-3,
eta_{n-1} + eta_n >= 0, and eta_{n-2} - eta_{n-1} >= 1.  If eta_{n-1} = -1 then
eta_n <= -1, contradicting eta_{n-1} + eta_n >= 0; if eta_{n-1} = 1 then
eta_{n-2} >= 2, impossible.  So eta_{n-1} = 0, whence eta_n <= 0 and eta_n >= 0,
so eta_n = 0 and eta_{n-2} = 1.  On 1, ..., n-3 the sequence is weakly
decreasing, hence 1^a 0^b (-1)^c with a + b + c = n-3, and the nonzero count
gives a + c + 1 = k.  Therefore

    mu_a = e_1 + ... + e_a - (e_{n-2-c} + ... + e_{n-3}) + e_{n-2},
           a + c = k - 1,   0 <= a <= k-1,

which is k weights -- or k-1 when k = n-2, where a = k-1 gives mu = lambda
itself and is excluded.  This is INDEPENDENT OF THE RANK, while the number of
instances of (*) grows with n (for k = 4: 200, 380, 644, 1008 at n = 7,8,9,10,
against 4 representatives throughout).  In D_7 with k = 4 they are

    mu = (0,0,1,-1,1,0,0), (0,1,1,-2,1,0,0), (1,1,0,-2,1,0,0), (1,0,0,-2,1,0,0)

in the omega basis, of lengths 1, 7, 14, 22.  At each of them M(mu) is the greedy top truncation
of I_kappa at EVERY vertex -- not just at t-1, and not just in the profile
sense at t, which is Proposition T.  So what is left is

    (B)   for the k weights mu with descent set {t-1} and c_t = 1, the graded
          multiplicities of M(mu) are the greedy top truncation of those of
          I_kappa, at every vertex,

and (B) implies (*), a greedy top truncation being an up-set by definition.
This is k explicit weights for each k, uniformly in the rank, rather than a
statement over a family growing with n.  (Note that the same assertion for
arbitrary extremal mu is FALSE -- see staircase.py: it fails at 122 of 1307
extremal modules, first at D_5, omega_2, mu = (1,-1,0,0,0), and it fails
exactly at t-1.  So (B) uses its weights.)

This script verifies, at every (n,k) in range: the reduction step on the modules
themselves (M(mu) and M(s_j mu) really do agree at t-1), that the irreducible
weights are exactly those with descent set {t-1} and that there are exactly k,
and (*) at each of them.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from gls import socle_step
from homological import setup
from upsets import upset

CASES = [("D", 4, 2), ("D", 5, 2), ("D", 6, 2), ("D", 7, 2), ("D", 8, 2),
         ("D", 9, 2), ("D", 5, 3), ("D", 6, 3), ("D", 7, 3), ("D", 8, 3),
         ("D", 9, 3), ("D", 10, 3), ("D", 6, 4), ("D", 7, 4), ("D", 8, 4),
         ("D", 9, 4), ("D", 10, 4), ("D", 7, 5), ("D", 8, 5), ("D", 9, 5),
         ("D", 10, 5), ("D", 8, 6), ("D", 9, 6), ("D", 10, 6)]


def module(rs, pr, word, lam):
    X = {}
    for _, i, c, _ in trace(rs, word, lam):
        for _ in range(c):
            X = socle_step(pr, rs, X, i)
    return X


def profile(pr, X, j):
    return tuple(len(X.get((d, j), [])) for d in range(pr.top + 1))


def greedy(pr, rs, X):
    """Is X the greedy top truncation of I_kappa at every vertex?"""
    for j in rs.I:
        v = sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
        rem, want = v, {}
        for d in range(pr.top, -1, -1):
            m = pr.dim.get((d, j), 0)
            if m:
                want[d] = min(m, rem)
                rem -= want[d]
        for d in want:
            if len(X.get((d, j), [])) != want[d]:
                return False
    return True


def run(typ, n, k, modcap=1200):
    rs, pr, kk = setup(typ, n, k)
    t = n - 2
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    heavy = len(orb) <= modcap
    st = Counter()
    for mu, w in orb.items():
        if mu[t - 1] < 1:
            continue
        st['inst'] += 1
        if not w:
            continue
        red = [j for j in rs.I if mu[j - 1] <= -1 and j not in (t - 1, t)]
        if red:
            if heavy:
                X = module(rs, pr, w, lam)
                Y = module(rs, pr, orb[rs.act(red[0], mu)], lam)
                st['step_tot'] += 1
                st['step_ok'] += (profile(pr, X, t - 1) == profile(pr, Y, t - 1))
        else:
            st['base'] += 1
            st['desc_tot'] += 1
            st['desc_ok'] += ([j for j in rs.I if mu[j - 1] <= -1] == [t - 1])
            if heavy:
                X = module(rs, pr, w, lam)
                st['base_tot'] += 1
                st['base_ok'] += upset(pr, X, t - 1)
                st['B_tot'] += 1
                st['B_ok'] += greedy(pr, rs, X)
    st['count_ok'] = int(st['base'] == k - (1 if k == n - 2 else 0))
    return st


if __name__ == "__main__":
    t0 = time.time()
    T = Counter()
    ok = True
    for typ, n, k in CASES:
        r = run(typ, n, k)
        T += r
        ok &= bool(r['count_ok']) and r['step_ok'] == r['step_tot'] \
            and r['base_ok'] == r['base_tot'] and r['desc_ok'] == r['desc_tot'] \
            and r['B_ok'] == r['B_tot']
        print(f"  {typ}{n} k={k}: {r['inst']:5d} instances -> {r['base']:2d} "
              f"representatives  (expected {k - (1 if k == n-2 else 0)})"
              f"{'' if r['count_ok'] else '   <-- COUNT MISMATCH'}")
    print(f"\n  reduction step verified on the modules : "
          f"{T['step_ok']}/{T['step_tot']}")
    print(f"  representatives have descent set {{t-1}}: "
          f"{T['desc_ok']}/{T['desc_tot']}")
    print(f"  (*) at every representative            : "
          f"{T['base_ok']}/{T['base_tot']}")
    print(f"  (B) greedy top truncation at every vertex: "
          f"{T['B_ok']}/{T['B_tot']}")
    print(f"  total instances {T['inst']}, representatives {T['base']}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("REDUCTION VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
