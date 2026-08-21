"""
THE STRAND DECOMPOSITION (W) -- proved, and with it (T0), (T1), (T2).

tau, the diagram automorphism exchanging the two leaves, fixes I_kappa when
kappa is a chain vertex.  At every CHAIN vertex it therefore splits the graded
pieces into eigenspaces, and the claim is that this splitting is as clean as it
could possibly be:

  (W1)  each eigenspace has multiplicity at most 1 in every chain slot, so a
        multiplicity-2 slot is exactly one + line and one - line;
  (W2)  the two eigen-multiplicities are given in closed form.  With t = n-2 and
        top degree 2t,

            a^+_d(j) = 1  <=>  |j-k| <= d <= j+k-2  and  d = j+k (mod 2),
            a^-_d(j) = a^+_{2t-d}(j) ,

        so the + strand is the Dirichlet half-line wave of the generator and the
        - strand is its reflection in the top degree;
  (W3)  each strand occupies a contiguous run of degrees (step 2) at every
        chain vertex: the + strand the low degrees, the - strand the high ones,
        overlapping exactly in the multiplicity-2 window;
  (W4)  every arrow restricts to an ISOMORPHISM between consecutive nonzero
        slots of the same strand.

So I_kappa is two multiplicity-free "waves" travelling through the chain, and
the arrows transport each wave.  Everything the conjecture still needs follows:
the rank of an arrow is the number of strands common to its two slots, which is
(T0); a rank-1 arrow into a multiplicity-2 slot carries whichever strand its
source has, and (T1) is the statement that this is the + one; and (T2) says the
+ strand has already died at a 1-dimensional chain slot one degree up.  Given
(W), all three are combinatorics of the two intervals, and thresh.py checks them
in that form over D_5--D_60 for every k.

WHAT IS PROVED.  All four.  (W1) and (W3) are immediate from (W2); (W2) follows
from the tau-refinement of Lemma R; and (W4) follows from the preprojective
relation at a chain vertex,

    a_{(j-1)->j} a_{j->(j-1)}  =  a_{(j+1)->j} a_{j->(j+1)}   on  I_j[d],

which on a thin strand turns the death of one arrow into the death of another
one vertex over, always moving AWAY from the arrow that died.  The step that
makes the induction go is the interval identity

    min( hi^eps_{j-1}, hi^eps_{j+1} )  =  hi^eps_j - 1,

immediate for the + strand, where hi^+_j = j+k-2, and for the - strand because
max(|j-1-k|, |j+1-k|) = |j-k|+1.  A death towards j+1 then walks down to vertex
1, which has a single chain neighbour, contradicting the tau-refined Lemma R; a
death towards j-1 walks up to t, where the relation carries the two leaf terms
and one needs

  (W4l)  the folded-leaf arrow  J[d]^eps -> I_t[d+1]^eps,  J = I_{n-1} (+) I_n,
         is nonzero whenever both sides are.

(W4l) is proved in turn: the leaves are multiplicity-free, so the image is
(1 + eps*tau)z with z = a_{(n-1)->t}u nonzero, and the relation at the leaf
gives a_{t->(n-1)}z = 0.  If I_t[d+1] is a line this is already the claim; if it
is 2-dimensional and z were an eigenvector, applying tau kills z with the other
leaf arrow too, so a_{t->(t-1)}z is nonzero and its outward arrow dies with the
target still in the strand -- contradicting the outward case, proved outright.

This script checks (W1)--(W4) and the closed form against the modules over
D_5--D_17 and every kappa, the arithmetic steps used in the proofs, and
(T0)--(T2) directly from the closed form over D_5--D_100.
"""
import sys, time
from collections import Counter
from homological import setup
from qvar import rref, P
from canon import build_tau, apply_tau
from thresh import images

CASES = [("D", n, k) for n in range(5, 18) for k in range(2, n - 1)]


def eigendim(pr, rs, T, d, j, sign):
    m = pr.dim[(d, j)]
    rows = []
    for b in range(m):
        e = [1 if x == b else 0 for x in range(m)]
        w = apply_tau(T, d, j, e)
        rows.append([(w[x] - sign * e[x]) % P for x in range(m)])
    A = [[rows[b][x] for b in range(m)] for x in range(m)]
    R, piv = rref(A, m)
    return m - len(piv)


def eigenbasis(pr, rs, T, d, j, sign):
    m = pr.dim[(d, j)]
    rows = []
    for b in range(m):
        e = [1 if x == b else 0 for x in range(m)]
        w = apply_tau(T, d, j, e)
        rows.append([(w[x] - sign * e[x]) % P for x in range(m)])
    A = [[rows[b][x] for b in range(m)] for x in range(m)]
    R, piv = rref(A, m)
    out = []
    for f in [c for c in range(m) if c not in piv]:
        x = [0] * m
        x[f] = 1
        for r, c in zip(R, piv):
            x[c] = (-r[f]) % P
        out.append(x)
    return out


def wave(n, k, j, d, eps, im=None):
    """a^eps_d(j) in {0,1}, in closed form."""
    if eps < 0:
        d = 2 * (n - 2) - d
    return 1 if (abs(j - k) <= d < j + k and (d - j - k) % 2 == 0) else 0


def arithmetic(nmax=101):
    """Nesting (hence (T0)), (T1) and (T2), straight from the closed form."""
    st = Counter()
    for n in range(5, nmax):
        t = n - 2
        for k in range(2, t + 1):
            for j in range(1, t + 1):
                for i in (j - 1, j + 1):
                    if not 1 <= i <= t:
                        continue
                    for d in range(0, 2 * t):
                        a1, a2 = wave(n, k, j, d, 1), wave(n, k, j, d, -1)
                        b1, b2 = wave(n, k, i, d + 1, 1), wave(n, k, i, d + 1, -1)
                        if not (a1 + a2) or not (b1 + b2):
                            continue
                        st['nest_tot'] += 1
                        st['nest_ok'] += ((a1 <= b1 and a2 <= b2)
                                          or (b1 <= a1 and b2 <= a2))
                        if a1 + a2 == 1 and b1 + b2 == 2:
                            st['T1_tot'] += 1
                            st['T1_ok'] += (a1 == 1)
                        if a1 + a2 == 2 and b1 + b2 == 1:
                            st['T2_tot'] += 1
                            st['T2_ok'] += (b2 == 1)
    return st


def run(typ, n, k):
    rs, pr, kk = setup(typ, n, k)
    t = n - 2
    T = build_tau(pr, rs, n, kk)
    im = None
    st = Counter()
    sup = {}
    for j in range(1, t + 1):
        for d in range(pr.top + 1):
            m = pr.dim.get((d, j), 0)
            wp = wave(n, k, j, d, 1, im)
            wm = wave(n, k, j, d, -1, im)
            st['W2_tot'] += 1
            st['W2_ok'] += (wp + wm == m and wp in (0, 1) and wm in (0, 1))
            if not m:
                continue
            ep, em = eigendim(pr, rs, T, d, j, 1), eigendim(pr, rs, T, d, j, -1)
            st['W1_tot'] += 1
            st['W1_ok'] += (ep <= 1 and em <= 1 and ep + em == m)
            st['W2e_tot'] += 1
            st['W2e_ok'] += (ep == wp and em == wm)
            for s, e in ((1, ep), (-1, em)):
                if e:
                    sup.setdefault((j, s), []).append(d)
    for (j, s), ds in sup.items():
        st['W3_tot'] += 1
        st['W3_ok'] += (ds == list(range(ds[0], ds[-1] + 1, 2)))
    for j in range(1, t + 1):
        for i in rs.adj[j]:
            if i > t:
                continue
            for d in range(pr.top):
                if not pr.dim.get((d, j)) or not pr.dim.get((d + 1, i)):
                    continue
                for s in (1, -1):
                    Bj = eigenbasis(pr, rs, T, d, j, s)
                    Bi = eigenbasis(pr, rs, T, d + 1, i, s)
                    if not Bj or not Bi:
                        continue
                    st['W4_tot'] += 1
                    st['W4_ok'] += any(y % P for y in pr.act(d, j, Bj[0], i))
    # duality: the rank identity predicted by D(I_kappa) = I_{kappa^o}
    def rk(d, i, j):
        mi, mj = pr.dim.get((d, i), 0), pr.dim.get((d + 1, j), 0)
        if not mi or not mj or (d, i, j) not in pr.arrow:
            return 0
        IM = [pr.act(d, i, [1 if x == c else 0 for x in range(mi)], j)
              for c in range(mi)]
        return len(rref([r[:] for r in IM], mj)[0])
    for j in rs.I:
        for i in rs.adj[j]:
            for d in range(pr.top):
                st['dual_tot'] += 1
                st['dual_ok'] += (rk(d, j, i) == rk(pr.top - d - 1, i, j))
    # the interval identity behind Lemma "a death propagates", both strands
    for j in range(2, t):
        for e in (1, -1):
            hi = lambda v: max([d for d in range(2 * t + 1)
                                if wave(n, k, v, d, e)], default=None)
            h0, hm, hp = hi(j), hi(j - 1), hi(j + 1)
            if None in (h0, hm, hp):
                continue
            st['hi_tot'] += 1
            st['hi_ok'] += (min(hm, hp) == h0 - 1)
    # the step used in (W4l): 2-dimensional slot at t forces (d+3,t) in strand
    for d in range(pr.top):
        if not wave(n, k, t + 1, d, 1):
            continue
        if not (wave(n, k, t, d + 1, 1) and wave(n, k, t, d + 1, -1)):
            continue
        for e in (1, -1):
            if not wave(n, k, t - 1, d + 2, e):
                continue
            st['c2_tot'] += 1
            st['c2_ok'] += bool(wave(n, k, t, d + 3, e))
    # (W4l): the folded-leaf arrow on each eigen-part
    for d in range(pr.top):
        b = pr.dim.get((d, n - 1), 0)
        m = pr.dim.get((d + 1, t), 0)
        if not b or not m:
            continue
        for s in (1, -1):
            if not eigenbasis(pr, rs, T, d + 1, t, s):
                continue
            for bi in range(b):
                u = [1 if x == bi else 0 for x in range(b)]
                i1 = pr.act(d, n - 1, u, t) if (d, n - 1, t) in pr.arrow else []
                i2 = (pr.act(d, n, apply_tau(T, d, n - 1, u), t)
                      if (d, n, t) in pr.arrow else [])
                img = [((i1[x] if x < len(i1) else 0)
                        + s * (i2[x] if x < len(i2) else 0)) % P
                       for x in range(m)]
                st['W4l_tot'] += 1
                st['W4l_ok'] += any(y % P for y in img)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        S += run(typ, n, k)
    S += arithmetic()
    rows = [("(W1) each strand multiplicity-free", 'W1_ok', 'W1_tot'),
            ("(W2) closed form = profile", 'W2_ok', 'W2_tot'),
            ("(W2) closed form = tau-eigendims", 'W2e_ok', 'W2e_tot'),
            ("(W3) each strand a contiguous interval", 'W3_ok', 'W3_tot'),
            ("(W4) arrows are isomorphisms on a strand", 'W4_ok', 'W4_tot'),
            ("(W4l) the folded-leaf arrow, each eigen-part", 'W4l_ok', 'W4l_tot'),
            ("duality: rank(j->i,d) = rank(i->j,2t-d-1)", 'dual_ok', 'dual_tot'),
            ("interval identity min(hi_{j-1},hi_{j+1}) = hi_j - 1",
             'hi_ok', 'hi_tot'),
            ("(W4l): the step (d+3,t) carries eps", 'c2_ok', 'c2_tot'),
            ("strand sets nested => (T0), to rank 100", 'nest_ok', 'nest_tot'),
            ("(T1) from the closed form, to rank 100", 'T1_ok', 'T1_tot'),
            ("(T2) from the closed form, to rank 100", 'T2_ok', 'T2_tot')]
    ok = True
    for lbl, a, b in rows:
        ok &= S[a] == S[b]
        print(f"  {lbl:42} {S[a]}/{S[b]}   "
              f"{'OK' if S[a] == S[b] else 'FAIL'}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("STRAND DECOMPOSITION VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
