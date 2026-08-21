"""
THE LAST STEP OF (B), AS ARITHMETIC -- AND TESTED AT SCALE.

canon.py reduces (B) to: the canonical top truncation T(v_a) is closed under
the arrows.  Three facts make every arrow of I_kappa transparent, each verified
there without exception:

  (T0)  every arrow between graded slots has rank min(dim source, dim target)
        -- it is injective or surjective, whichever is possible (450/450);
  (T1)  a rank-1 arrow into a multiplicity-2 CHAIN slot has image the
        tau-eigenline (22/22);
  (T2)  the tau-eigenline of a multiplicity-2 chain slot is killed by exactly
        the arrows into 1-dimensional chain slots (102/102).

So whether a closure obligation is satisfied depends only on the DIMENSIONS of
the two slots, on whether the target is a leaf, and on the thresholds.  That
makes the remaining content of (B) pure arithmetic, and this script tests it far
beyond the range where the modules themselves can be built.

CLOSED FORMS.  For D_n, lambda = omega_k, k <= n-2, t = n-2, the injective is
Pi e_k and its top degree is 2n-4.  At a vertex j the occupied degrees are
    dist(k,j), dist(k,j)+2, ..., (2n-4) - dist(k,j),
with multiplicity 2 in the middle max(0, j+k-n+1) of them for a chain vertex
j <= t, and multiplicity 1 throughout at the two leaves (Lemma P).  For the
representative mu_a of Proposition (reps), in Bourbaki coordinates
    x = lambda - mu_a,   v_j = x_1 + ... + x_j  (j <= n-2),
    v_{n-1} = v_n = v_{n-2}/2,
since x vanishes in the last two coordinates.

Both closed forms are checked against the actual module in small rank before
being used, and then the obligation classification is run over a large range.
"""
import sys, time
from collections import Counter

CHAIN, LEAF = 'chain', 'leaf'


def images(n, k, span=6):
    """The images of k under the infinite dihedral group of Lemma P, signed.

    Dirichlet at 0 (odd reflection, sign flip) and Neumann at L = n-1 (even
    reflection, sign kept), starting from k with sign +1.
    """
    L = n - 1
    out = {}
    frontier = [(k, 1)]
    out[k] = 1
    for _ in range(span):
        new = []
        for x, s in frontier:
            for y, t in ((-x, -s), (2 * L - x, s)):
                if abs(y) <= 4 * L * span and y not in out:
                    out[y] = t
                    new.append((y, t))
        frontier = new
    return out


def profile(n, k, j, _c={}, _p={}):
    """[(degree, multiplicity)] of I_kappa at vertex j, by the method of images.

    a_d(j) = sum over signed images x of  sign(x) * [ |j-x| <= d, j = x+d (2) ].
    Multiplicities can be 0 inside the range, not only 1 or 2, so the occupied
    degrees are read off rather than assumed to be a full progression.
    """
    if (n, k, j) in _p:
        return _p[(n, k, j)]
    key = (n, k)
    if key not in _c:
        _c[key] = images(n, k)
    im = _c[key]
    top = 2 * n - 4
    t = n - 2
    # the folding of Lemma P: the two leaves share the value at index t+1 = n-1,
    # each carrying half of it
    idx, half = (j, False) if j <= t else (n - 1, True)
    out = []
    for d in range(top + 1):
        v = 0
        for x, s in im.items():
            if abs(idx - x) <= d and (idx - x - d) % 2 == 0:
                v += s
        if half:
            v //= 2
        if v:
            out.append((d, v))
    _p[(n, k, j)] = out
    return out


def wave(n, k, j, d, eps, _c={}):
    """a^eps_d(j) in {0,1}: the eps-strand of I_kappa at the chain slot (d,j).

    The signed images of Lemma P split into two classes modulo 4L, L = n-1:
    those congruent to +-k -- the strand of the generator -- and the rest, the
    strand of its mirror image at L.  Each class contributes 0 or 1, and the two
    add up to the multiplicity.  The classes are the tau-eigenspaces: the +k
    class is the tau-invariant one.
    """
    key = (n, k)
    if key not in _c:
        _c[key] = images(n, k)
    L = n - 1
    v = 0
    for x, s in _c[key].items():
        if abs(j - x) <= d and (j - x - d) % 2 == 0:
            near = ((x - k) % (4 * L) == 0 or (x + k) % (4 * L) == 0)
            if near == (eps > 0):
                v += s
    return v


def strand_counts(n, k, j, st):
    """|S^eps(j)| = min(j,k) and the overlap = max(0, j+k-t-1)."""
    t = n - 2
    sp = sum(1 for d in range(2 * t + 1) if wave(n, k, j, d, 1))
    sm = sum(1 for d in range(2 * t + 1) if wave(n, k, j, d, -1))
    ov = sum(1 for d in range(2 * t + 1)
             if wave(n, k, j, d, 1) and wave(n, k, j, d, -1))
    st['cnt_tot'] += 1
    st['cnt_ok'] += (sp == min(j, k) and sm == min(j, k)
                     and ov == max(0, j + k - t - 1))


def thresholds(n, k, a):
    """The greedy threshold D_j: the highest degree not taken in full."""
    v = vvec(n, k, a)
    D = {}
    for j in range(1, n - 1):
        rem, thr = v[j], -1
        for d, m in reversed(profile(n, k, j)):
            if rem >= m:
                rem -= m
                continue
            thr = d
            break
        D[j] = thr
    return D


def occupation(n, k, a, j):
    """(s^+_j, s^-_j): how many slots of each strand the truncation takes.

    With M = min(j,k) slots in each strand and o = max(0, j+k-t-1) of them
    shared, greedy from the top fills the M-o slots that carry only the minus
    strand, then the o shared ones (two dimensions each, the odd one taking the
    tau-line, which is the plus strand), then the plus-only slots.
    """
    t = n - 2
    M, o, v = min(j, k), max(0, j + k - t - 1), vvec(n, k, a)[j]
    if v <= M - o:
        return 0, v
    r = v - (M - o)
    if r <= 2 * o:
        return (r + 1) // 2, (M - o) + r // 2
    return v - M, M


def increments(n, k, a, st):
    """The occupation numbers move by at most one along the chain.

    Writing u_j = v_j - min(j,k) and o_j = max(0, j+k-t-1), the closed form is
    s^+_j = max(0, u_j, ceil((u_j+o_j)/2)), and the proof of the increment bounds
    rests on four facts, all checked here:

      (P1)  u_{j+1} - u_j is -1, 0 or 1;
      (P2)  it is -1 only for j+1 <= a, where s^+ vanishes at both ends, or for
            j+1 = t, where o steps up so that u+o is unchanged;
      (P3)  u_{t-1} <= o_{t-1} + 1, so at j+1 = t the term ceil((u+o)/2) is not
            strictly beaten by u and s^+ cannot drop;
      (P5)  where the dimension vector jumps by 2, s^+ rises by exactly 1 --
            either u already attains the maximum, or o steps up and carries
            ceil((u+o)/2) with it.  These two account for every instance; the
            parity case that would break the argument never occurs.
    """
    t = n - 2
    v = vvec(n, k, a)
    u = lambda j: v[j] - min(j, k)
    o = lambda j: max(0, j + k - t - 1)
    st['P3_tot'] += 1
    st['P3_ok'] += (u(t - 1) <= o(t - 1) + 1)
    for j in range(1, t):
        sp0, sm0 = occupation(n, k, a, j)
        sp1, sm1 = occupation(n, k, a, j + 1)
        st['inc_tot'] += 1
        st['inc_ok'] += (0 <= sp1 - sp0 <= 1)
        st['incm_tot'] += 1
        st['incm_ok'] += (abs(sm1 - sm0) <= 1)
        du = u(j + 1) - u(j)
        st['P1_tot'] += 1
        st['P1_ok'] += (-1 <= du <= 1)
        if du == -1:
            st['P2_tot'] += 1
            st['P2_ok'] += (j + 1 <= a or j + 1 == t)
        if v[j + 1] - v[j] == 2:
            st['P5_tot'] += 1
            st['P5_ok'] += (sp1 - sp0 == 1)
            C = -((-(u(j) + o(j))) // 2)
            if u(j) >= C:
                st['P5_byu'] += 1
            elif o(j + 1) > o(j):
                st['P5_ostep'] += 1
            else:
                st['P5_gap'] += 1


def upclosed(n, k, a, st):
    """The truncation, split by strand, is up-closed in each strand; and the
    thresholds move by at most one along the chain away from the node."""
    t = n - 2
    v = vvec(n, k, a)
    U = {1: set(), -1: set()}
    for j in range(1, t + 1):
        rem = v[j]
        cnt = {1: 0, -1: 0}
        for d, m in reversed(profile(n, k, j)):
            if rem <= 0:
                break
            take = min(m, rem)
            rem -= take
            if take == m:
                for e in (1, -1):
                    if wave(n, k, j, d, e):
                        U[e].add((d, j))
                        cnt[e] += 1
            else:
                U[1].add((d, j))          # the tau-line is the + line
                cnt[1] += 1
        st['occ_tot'] += 1
        st['occ_ok'] += ((cnt[1], cnt[-1]) == occupation(n, k, a, j))
    for e in (1, -1):
        for (d, j) in U[e]:
            for i in (j - 1, j + 1):
                if not 1 <= i <= t or not wave(n, k, i, d + 1, e):
                    continue
                st['up_tot'] += 1
                st['up_ok'] += ((d + 1, i) in U[e])
    D = thresholds(n, k, a)
    for j in range(1, t - 1):
        st['adj_tot'] += 1
        st['adj_ok'] += (abs(D[j + 1] - D[j]) <= 1)


def leaf_boundary(n, k, a, st):
    """The obligations that involve a leaf always land in a full slot.

    At the trivalent node v_t = 2k-2a-2 and each leaf carries v = k-a-1, so the
    truncation at t has lowest taken degree, and threshold, equal to t-k+2a+2,
    while the leaves have lowest taken degree t-k+2a+3 -- exactly one more.  So
    every taken slot at t maps into a taken leaf slot, and every taken leaf slot
    maps into the slot of degree (t-k+2a+2)+2 at t, which is above the threshold
    and hence full.
    """
    t = n - 2
    v = vvec(n, k, a)
    Tr = truncation(n, k, v)
    if Tr is None:
        return
    st['Q1_tot'] += 1
    st['Q1_ok'] += (v[t] == 2 * k - 2 * a - 2 and v[n - 1] == k - a - 1)
    prt, prl = profile(n, k, t), profile(n, k, n - 1)
    tk_t = [d for d, _ in prt if Tr[t].get(d, 0) > 0]
    tk_l = [d for d, _ in prl if Tr[n - 1].get(d, 0) > 0]
    if v[t] >= 2:
        D = max([d for d, m in prt if Tr[t].get(d, 0) < m], default=-1)
        st['Q2_tot'] += 1
        st['Q2_ok'] += (min(tk_t) == t - k + 2 * a + 2 == D)
    if v[n - 1] >= 1:
        st['Q3_tot'] += 1
        st['Q3_ok'] += (min(tk_l) == t - k + 2 * a + 3)
    dt, dl = dict(prt), dict(prl)
    for d in tk_t:
        for L in (n - 1, n):
            if d + 1 in dl:
                st['Q4_tot'] += 1
                st['Q4_ok'] += (Tr[L].get(d + 1, 0) == dl[d + 1])
    for d in tk_l:
        if d + 1 in dt:
            st['Q4_tot'] += 1
            st['Q4_ok'] += (Tr[t].get(d + 1, 0) == dt[d + 1])


def local_facts(n, k, st):
    """(T0), (T1), (T2) as consequences of the strand decomposition (W)."""
    t = n - 2
    for j in range(1, t + 1):
        dj = dict(profile(n, k, j))
        for d, mj in dj.items():
            for i in neighbours(n, j):
                mi = dim_at(n, k, i, d + 1)
                if not mi:
                    continue
                if i > t:                      # a leaf: tau exchanges the two
                    continue
                r = sum(1 for e in (1, -1)
                        if wave(n, k, j, d, e) and wave(n, k, i, d + 1, e))
                st['T0_tot'] += 1
                st['T0_ok'] += (r == min(mj, mi))
                if r == 1 and mi == 2:
                    st['T1_tot'] += 1
                    st['T1_ok'] += bool(wave(n, k, j, d, 1)
                                        and wave(n, k, i, d + 1, 1))
                if mj == 2 and mi == 1:
                    st['T2_tot'] += 1
                    st['T2_ok'] += (wave(n, k, i, d + 1, 1) == 0)


def dim_at(n, k, j, d):
    return dict(profile(n, k, j)).get(d, 0)


def vvec(n, k, a):
    """dim V(mu_a) as a root-coordinate vector, indices 1..n."""
    c = k - 1 - a
    b = n - 2 - k
    x = [0] * (n + 1)
    for p in range(1, n + 1):
        lam = 1 if p <= k else 0
        if p <= a:
            m = 1
        elif p <= a + b:
            m = 0
        elif p <= n - 3:
            m = -1
        elif p == n - 2:
            m = 1
        else:
            m = 0
        x[p] = lam - m
    v = [0] * (n + 1)
    s = 0
    for p in range(1, n - 1):
        s += x[p]
        v[p] = s
    v[n - 1] = v[n] = v[n - 2] // 2
    return v


def truncation(n, k, v):
    """{j: {degree: taken}} for the greedy top truncation with totals v."""
    out = {}
    for j in range(1, n + 1):
        pr = profile(n, k, j)
        rem, row = v[j], {}
        for d, m in reversed(pr):
            take = min(m, rem)
            rem -= take
            row[d] = take
        if rem:
            return None
        out[j] = row
    return out


def neighbours(n, j):
    t = n - 2
    if j <= t:
        return [i for i in (j - 1, j + 1) if 1 <= i <= t] + \
               ([n - 1, n] if j == t else [])
    return [t]


def classify(n, k, a):
    """Classify every closure obligation of T(v_a).  Returns a Counter."""
    v = vvec(n, k, a)
    T = truncation(n, k, v)
    st = Counter()
    if T is None:
        st['BAD_dimension'] += 1
        return st
    dim = {j: dict(profile(n, k, j)) for j in range(1, n + 1)}
    for j in range(1, n + 1):
        for d, take in T[j].items():
            if take == 0:
                continue
            src_full = (take == dim[j][d])
            for i in neighbours(n, j):
                m = dim[i].get(d + 1, 0)
                if not m:
                    continue
                got = T[i].get(d + 1, 0)
                tgt = 'full' if got == m else ('partial' if got else 'empty')
                leaf = i > n - 2
                rank = min(dim[j][d], m)          # (T0)
                if src_full:
                    if tgt == 'full':
                        st['full->full'] += 1
                    elif tgt == 'partial':
                        # target has m = 2, got = 1: need the image inside the
                        # eigenline, so the arrow must have rank 1, and then
                        # (T1) identifies the image.
                        if rank == 1 and not leaf:
                            st['full->partial (T1)'] += 1
                        else:
                            st['BAD full->partial rank2'] += 1
                    else:
                        st['BAD full->empty'] += 1
                else:
                    # the source is the tau-eigenline of a 2-slot
                    if tgt == 'full':
                        st['line->full'] += 1
                    elif tgt == 'partial':
                        st['line->partial (equivariance)'] += 1
                    else:
                        if m == 1 and not leaf:
                            st['line->empty (T2)'] += 1
                        else:
                            st['BAD line->empty'] += 1
    return st


def check_closed_forms():
    """Validate profile() and vvec() against the actual modules."""
    from homological import setup
    from rootsys import orbit_with_words
    import reduce as R
    st = Counter()
    for n, k in [(5, 2), (6, 2), (7, 2), (6, 3), (7, 3), (8, 3), (7, 4),
                 (8, 4), (8, 5)]:
        rs, pr, kk = setup("D", n, k)
        st['kk_tot'] += 1
        st['kk_ok'] += (kk == k)
        for j in rs.I:
            got = [(d, pr.dim[(d, j)]) for d in range(pr.top + 1)
                   if pr.dim.get((d, j))]
            st['prof_tot'] += 1
            st['prof_ok'] += (got == profile(n, k, j))
        t = n - 2
        lam = rs.fundamental(k)
        vs = set()
        for mu, w in orbit_with_words(rs, lam).items():
            if mu[t - 1] < 1 or not w:
                continue
            if any(mu[j - 1] <= -1 and j not in (t - 1, t) for j in rs.I):
                continue
            X = R.module(rs, pr, w, lam)
            vs.add(tuple(sum(len(X.get((d, j), []))
                             for d in range(pr.top + 1)) for j in rs.I))
        pred = set()
        for a in range(k):
            v = vvec(n, k, a)
            if any(x < 0 for x in v[1:]):
                continue
            pred.add(tuple(v[1:]))
        pred.discard(tuple([0] * n))
        st['v_tot'] += 1
        st['v_ok'] += (vs == pred)
    return st


if __name__ == "__main__":
    t0 = time.time()
    C = check_closed_forms()
    print(f"  kappa' = k                {C['kk_ok']}/{C['kk_tot']}")
    print(f"  Lemma P profile matches   {C['prof_ok']}/{C['prof_tot']}")
    print(f"  v_a closed form matches   {C['v_ok']}/{C['v_tot']}")
    S = Counter()
    N = 0
    for n in range(5, 46):
        for k in range(2, n - 1):
            for a in range(k):
                v = vvec(n, k, a)
                if any(x < 0 for x in v[1:]) or all(x == 0 for x in v[1:]):
                    continue
                S += classify(n, k, a)
                upclosed(n, k, a, S)
                increments(n, k, a, S)
                leaf_boundary(n, k, a, S)
                N += 1
            local_facts(n, k, S)
            for j in range(1, n - 1):
                strand_counts(n, k, j, S)
    bad = {a: b for a, b in S.items() if a.startswith('BAD')}
    print(f"\n  closure obligations over {N} representatives "
          f"(D_5..D_45, every k, every a):")
    for a, b in sorted(S.items()):
        if not a.startswith('BAD') and not a.startswith(('T0_','T1_','T2_','up_','adj_','cnt_','occ_','inc','P1_','P2_','P3_','P5','Q1_','Q2_','Q3_','Q4_')):
            print(f"      {a:32} {b}")
    ok = (not bad and C['kk_ok'] == C['kk_tot']
          and C['prof_ok'] == C['prof_tot'] and C['v_ok'] == C['v_tot']
          and S['T0_ok'] == S['T0_tot'] and S['T1_ok'] == S['T1_tot']
          and S['T2_ok'] == S['T2_tot'] and S['cnt_ok'] == S['cnt_tot']
          and S['up_ok'] == S['up_tot'] and S['adj_ok'] == S['adj_tot']
          and S['occ_ok'] == S['occ_tot'] and S['inc_ok'] == S['inc_tot']
          and S['incm_ok'] == S['incm_tot'] and S['P1_ok'] == S['P1_tot']
          and S['P2_ok'] == S['P2_tot'] and S['P3_ok'] == S['P3_tot']
          and S['P5_ok'] == S['P5_tot'] and S['P5_gap'] == 0
          and S['Q1_ok'] == S['Q1_tot'] and S['Q2_ok'] == S['Q2_tot']
          and S['Q3_ok'] == S['Q3_tot'] and S['Q4_ok'] == S['Q4_tot'])
    if bad:
        print("      UNRESOLVED:", bad)
    print(f"\n  (T0) rank = min(dim,dim)      {S['T0_ok']}/{S['T0_tot']}")
    print(f"  (T1) rank-1 image is the + strand  {S['T1_ok']}/{S['T1_tot']}")
    print(f"  (T2) + strand dies into 1-dim chain slots  "
          f"{S['T2_ok']}/{S['T2_tot']}")
    print(f"\n  |S^eps(j)| = min(j,k), overlap = max(0,j+k-t-1)  "
          f"{S['cnt_ok']}/{S['cnt_tot']}")
    print(f"  the truncation is up-closed in each strand      "
          f"{S['up_ok']}/{S['up_tot']}")
    print(f"  |D_j - D_(j+1)| <= 1 away from the node         "
          f"{S['adj_ok']}/{S['adj_tot']}")
    print(f"  the occupation formula reproduces the truncation "
          f"{S['occ_ok']}/{S['occ_tot']}")
    print(f"  0 <= s+_(j+1) - s+_j <= 1                       "
          f"{S['inc_ok']}/{S['inc_tot']}")
    print(f"  |s-_(j+1) - s-_j| <= 1                          "
          f"{S['incm_ok']}/{S['incm_tot']}")
    print(f"    (P1) u increments in {{-1,0,1}}                  "
          f"{S['P1_ok']}/{S['P1_tot']}")
    print(f"    (P2) a drop only at p <= a or at t             "
          f"{S['P2_ok']}/{S['P2_tot']}")
    print(f"    (P3) u_(t-1) <= o_(t-1) + 1                    "
          f"{S['P3_ok']}/{S['P3_tot']}")
    print(f"    (P5) a jump of 2 raises s+ by exactly 1        "
          f"{S['P5_ok']}/{S['P5_tot']}"
          f"   [u attains it {S['P5_byu']}, o steps up {S['P5_ostep']}, "
          f"neither {S['P5_gap']}]")
    print(f"    (Q1) v_t = 2k-2a-2, v_leaf = k-a-1             "
          f"{S['Q1_ok']}/{S['Q1_tot']}")
    print(f"    (Q2) delta_t = D_t = t-k+2a+2                  "
          f"{S['Q2_ok']}/{S['Q2_tot']}")
    print(f"    (Q3) the leaves start one degree higher        "
          f"{S['Q3_ok']}/{S['Q3_tot']}")
    print(f"    (Q4) every leaf obligation lands in a full slot "
          f"{S['Q4_ok']}/{S['Q4_tot']}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("ARITHMETIC VERIFIED -- every obligation resolved by (T0),(T1),(T2)"
          if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
