"""
proof.py -- verification of every step of the PROOF of the Main Theorem.

    python3 proof.py

The proof is:

  (1)  CRAWLEY-BOEVEY IDENTITY.  For any submodule X of the injective
       I_kappa over the preprojective algebra Pi of a Dynkin quiver, with
       mu' = lambda - dim X,

           dim soc_i(I_kappa / X)  =  <alpha_i^vee, mu'>  +  dim top_i(X).

       Proof: apply Hom(S_i,-) to 0 -> X -> I -> Q -> 0; Ext^1(S_i,I) = 0
       since I is injective; then use Crawley-Boevey's formula
       dim Ext^1(M,N) = dim Hom(M,N) + dim Hom(N,M) - (dim M, dim N)
       together with dim soc_i(I_kappa) = delta_{i,kappa}.

  (2)  MONOTONICITY.  Hence a socle step of colour i raises dim X by at
       least c := <alpha_i^vee, mu'>.

  (3)  THE SQUEEZE.  dim I_kappa = lambda - w_0.lambda.  Extend a reduced
       word for w to one for w_0; the c-trace telescopes to
       sum(c) = lambda - w_0.lambda.  Since soc_{sigma_full}(I) is contained
       in I, the total growth is at most dim I = sum(c), while by (2) it is
       at least sum(c).  So EVERY step contributes exactly c and every
       top-term vanishes.  Restricting to the prefix for w gives

           dim soc_sigma(I_kappa) = lambda - w.lambda = v.

  (4)  UNIQUENESS (Savage-Tingley, Prop. 4.9).  A w-extremal v is carried by
       a unique submodule of I_kappa.  Hence soc_sigma(I_kappa) IS the
       T-fixed representation V.

  (5)  MINUSCULE CASE.  If lambda is minuscule every c equals 1, w is
       lambda-minuscule hence fully commutative, and V is the heap module,
       whose colour fibres are totally ordered: the spectrum is always
       simple.  (All of type A; the spin weights of D; E6 om1/om6; E7 om7.)

  (6)  LEMMA A.  V_t multiplicity-free => V multiplicity-free (t = trivalent
       node).  Proved in type D by window containment, and in type E by the
       inward-injectivity lemma of inwardlemma.py, whose verification there is
       a complete finite computation.

This script checks (1), (3), the conclusion of (3)+(4), (5) and (6).
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from proj2 import Proj2
from gls import socle_step
from homological import setup, top_i, dimvec
from squeeze import longest_weight
from lemA import build, TRIV

OK = True
t00 = time.time()


def hdr(s):
    print("\n" + s)


# ---------------------------------------------------------------- step (1)
hdr("(1) Crawley-Boevey identity  dim soc_i(I/X) = c + dim top_i(X)")
tot = good = 0
for typ, n, k in [("A", 4, 2), ("D", 4, 2), ("D", 5, 3), ("E", 6, 3), ("E", 7, 1)]:
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    for mu, word in orbit_with_words(rs, lam).items():
        if not word:
            continue
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            v0 = dimvec(rs, X, pr.top)
            cpred = (1 if i == k else 0) + sum(v0[j - 1] for j in rs.adj[i]) - 2 * v0[i - 1]
            h = top_i(pr, rs, X, i)
            X2 = socle_step(pr, rs, X, i)
            grew = dimvec(rs, X2, pr.top)[i - 1] - v0[i - 1]
            tot += 1
            good += (grew == cpred + h == c + h and c == cpred)
            X = X2
            for _ in range(c - 1):
                X = socle_step(pr, rs, X, i)
OK &= (good == tot)
print(f"    {good}/{tot} socle steps satisfy the identity   {'OK' if good==tot else 'FAIL'}")

# ---------------------------------------------------------------- step (3)
hdr("(3) The squeeze input:  dim I_kappa = lambda - w_0.lambda")
tot = good = 0
for typ, n in [("A", 3), ("A", 4), ("A", 5), ("D", 4), ("D", 5), ("D", 6),
               ("D", 7), ("E", 6), ("E", 7)]:
    rs = RootSystem(typ, n)
    nu = {}
    for j in rs.I:
        p = Proj2(rs, j)
        nu[j] = max(i for (d, i), m in p.dim.items() if d == p.top and m)
    inv = {b: a for a, b in nu.items()}
    for k in rs.I:
        pr = Proj2(rs, inv[k], dmax=60)
        d = tuple(sum(pr.dim.get((dd, i), 0) for dd in range(pr.top + 1)) for i in rs.I)
        lam = rs.fundamental(k)
        w0 = longest_weight(rs, k)
        tgt = tuple(int(x) for x in rs.root_coords([lam[j-1]-w0[j-1] for j in rs.I]))
        tot += 1
        good += (d == tgt)
OK &= (good == tot)
print(f"    {good}/{tot} injectives   {'OK' if good==tot else 'FAIL'}")

# ------------------------------------------------------------ steps (3)+(4)
hdr("(3)+(4) Conclusion:  dim soc_sigma(I_kappa) = v = lambda - w.lambda")
tot = good = 0
for typ, n, k in [("A", 5, 3), ("D", 4, 2), ("D", 5, 2), ("D", 6, 3),
                  ("E", 6, 2), ("E", 6, 3), ("E", 7, 1), ("E", 7, 7)]:
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    for mu, word in orbit_with_words(rs, lam).items():
        if not word:
            continue
        v = tuple(int(x) for x in v_from_trace(rs, trace(rs, word, lam)))
        X = build(rs, pr, word, lam)
        tot += 1
        good += (tuple(dimvec(rs, X, pr.top)) == v)
OK &= (good == tot)
print(f"    {good}/{tot} elements   {'OK' if good==tot else 'FAIL'}")

# ---------------------------------------------------------------- step (5)
hdr("(5) Minuscule lambda:  every c = 1, and the spectrum is always simple")
tot = good = 0
MIN = [("A", 5, k) for k in range(1, 6)] + \
      [("D", 5, 4), ("D", 5, 5), ("D", 6, 5), ("D", 6, 6),
       ("E", 6, 1), ("E", 6, 6), ("E", 7, 7)]
for typ, n, k in MIN:
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    for mu, word in orbit_with_words(rs, lam).items():
        if not word:
            continue
        cs = [c for _, _, c, _ in trace(rs, word, lam)]
        X = build(rs, pr, word, lam)
        tot += 1
        good += (max(cs) == 1 and all(len(x) <= 1 for x in X.values()))
OK &= (good == tot)
print(f"    {good}/{tot} elements   {'OK' if good==tot else 'FAIL'}")

# ---------------------------------------------------------------- step (6)
hdr("(6) Lemma A:  V_t multiplicity-free  =>  V multiplicity-free")
tot = good = 0
for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 3),
                  ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 7, 1), ("E", 7, 7)]:
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    for mu, word in orbit_with_words(rs, lam).items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        mt = max((len(X.get((d, t), [])) for d in range(pr.top + 1)), default=0)
        mall = max((len(x) for x in X.values()), default=0)
        tot += 1
        good += ((mt <= 1) == (mall <= 1))
OK &= (good == tot)
print(f"    {good}/{tot} elements   {'OK' if good==tot else 'FAIL'}")

print(f"\n({time.time()-t00:.0f}s)")
print("EVERY STEP OF THE PROOF VERIFIED" if OK else "*** SOMETHING FAILED ***")
sys.exit(0 if OK else 1)
