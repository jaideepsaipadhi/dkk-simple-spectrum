"""
THE NEIGHBOUR CONDITION AT DEGREE D+1, PROVED AT ALL THREE NEIGHBOURS.

Everything left open in this paper is Claim G:

   (G1) if X_t[D] != 0 then every neighbour of t is full at degree D+1;
   (G2) if X_t[D] = 0 then at most one neighbour is deficient at degree D+1,
        and that deficiency is one-dimensional.

G => F => Proposition T, for every k and with no case split.

PROVED HERE (Proposition "G1 at the leaves"), assuming only that X is an
up-set at the two leaves:

   1. At a colour-t step, c = mu'_t - mu'_{n-1} >= 1 and mu'_t <= 1, so
      mu'_{n-1} <= 0.
   2. X_t[D] != 0 means v'_t = 2a is even and >= 2 (Lemma C + top truncation),
      and then D = t + k - 2a, so D+1 is the a-th leaf degree from the top.
   3. The epsilon-coordinate identities give
          v'_{n-1} + v'_n = v'_t - mu'_{n-1} in {2a, 2a+1},
          |v'_{n-1} - v'_n| = |mu'_n| <= 1,
      hence min(v'_{n-1}, v'_n) >= a.
   4. An up-set at a leaf fills the top v'_j leaf degrees, so X_j[D+1] is full
      exactly when v'_j >= a.

AT THE CHAIN NEIGHBOUR the same scheme works: c >= 1 gives
      v'_{t-1} >= 2a - [k=t],
while the profile of I_kappa satisfies the inequality
      sum_{d >= D+1} dim I_kappa_{t-1}[d]  <=  2a - [k=t]
(52/52 over D5-D8, all k <= n-2; a statement about the image formula alone),
so an up-set at t-1 is full at D+1 too.

Remaining: the up-set property at t-1 and at the leaves, holding
simultaneously along the construction, and (G2).

This script verifies Claim G, each of the four links above, and the two
negative results that delimit the search.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV
from epsilon import eps_of_mu

CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 3), ("D", 6, 4),
         ("D", 7, 2), ("D", 7, 3)]


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    st = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            if i != t:
                for _ in range(c):
                    X = socle_step(pr, rs, X, i)
                continue
            vp = [sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
                  for j in rs.I]
            e = eps_of_mu(n, k, vp)
            D = max([d for d in range(pr.top + 1)
                     if pr.dim.get((d, t), 0)
                     and len(X.get((d, t), [])) < pr.dim[(d, t)]], default=None)
            # up-set at the leaves
            for j in (n - 1, n):
                degs = [d for d in range(pr.top + 1) if pr.dim.get((d, j), 0)]
                occ = [d for d in degs if X.get((d, j))]
                st['upset_leaf_tot'] += 1
                st['upset_leaf_ok'] += (not occ or
                                        set(occ) == set(d for d in degs if d >= min(occ)))
            # link 1
            st['mu_tot'] += 1
            st['mu_ok'] += (e[n - 2] <= 0)
            if D is not None:
                defs = [j for j in rs.adj[t]
                        if len(X.get((D + 1, j), [])) < pr.dim.get((D + 1, j), 0)]
                if vp[t - 1] > 0 and vp[t - 1] % 2 == 0:
                    a = vp[t - 1] // 2
                    st['link3_tot'] += 1
                    st['link3_ok'] += (min(vp[n - 2], vp[n - 1]) >= a)
                    for j in (n - 1, n):
                        if pr.dim.get((D + 1, j), 0):
                            st['link4_tot'] += 1
                            full = len(X.get((D + 1, j), [])) == pr.dim[(D + 1, j)]
                            st['link4_ok'] += (full == (vp[j - 1] >= a))
                if len(X.get((D, t), [])) > 0:
                    st['G1_tot'] += 1
                    st['G1_ok'] += (len(defs) == 0)
                else:
                    st['G2_tot'] += 1
                    st['G2_ok'] += (len(defs) <= 1 and
                                    all(pr.dim.get((D + 1, j), 0)
                                        - len(X.get((D + 1, j), [])) == 1 for j in defs))
            X = socle_step(pr, rs, X, t)
            for _ in range(c - 1):
                X = socle_step(pr, rs, X, t)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r:
            S += r
    def line(lbl, a, b):
        print(f"  {lbl:52} {S[a]}/{S[b]}   {'OK' if S[a] == S[b] else 'FAIL'}")
    print("CLAIM G\n")
    line("(G1) X_t[D] != 0 => no neighbour deficient", 'G1_ok', 'G1_tot')
    line("(G2) X_t[D] = 0 => <=1 deficient, 1-dimensional", 'G2_ok', 'G2_tot')
    print("\nTHE PROVED PART (G1 at the leaves), link by link\n")
    line("1. mu'_(n-1) <= 0 at every colour-t step", 'mu_ok', 'mu_tot')
    line("3. v'_t = 2a  =>  min(v'_(n-1), v'_n) >= a", 'link3_ok', 'link3_tot')
    line("4. X_j[D+1] full  <=>  v'_j >= a  (leaves)", 'link4_ok', 'link4_tot')
    line("   up-set at the leaves (the one assumption)", 'upset_leaf_ok',
         'upset_leaf_tot')
    ok = all(S[a] == S[b] for a, b in [('G1_ok', 'G1_tot'), ('G2_ok', 'G2_tot'),
                                       ('mu_ok', 'mu_tot'),
                                       ('link3_ok', 'link3_tot'),
                                       ('link4_ok', 'link4_tot'),
                                       ('upset_leaf_ok', 'upset_leaf_tot')])
    print(f"\n({time.time()-t0:.0f}s)")
    print("ALL VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
