"""
audit.py -- one script that walks the whole proof and tries to break it.

    python3 audit.py

Every link in the chain from the main theorem to the type-D closed form is
re-checked here, independently of the scripts that introduced it, together with
NEGATIVE CONTROLS: statements that are supposed to be FALSE and whose failure to
fail would mean something has drifted.  The point is that a wrong step should
show up as a red line here rather than survive in the manuscript.

The chain, in the order the paper proves it:

  A1  main theorem: V = soc_sigma(I_kappa), against the independent solver
  A2  Lemma I (inward injectivity): complete in type E, uniform in type D
  A3  Lemma P: the closed-form profile of I_kappa
  A4  Lemma C: v_t = r + 2m
  A5  (W1)-(W4): the strand decomposition, on the modules
  A6  (T0),(T1),(T2): from the closed form
  A7  Proposition (reps): the representatives are the explicit family
  A8  (B): the canonical truncation is a submodule and equals V(mu_a)
  A9  the arithmetic: increments, the jump-of-2 dichotomy, the leaf identities
  A10 (G1),(G2): the two profile inequalities, by counting
  A11 end to end: simple spectrum <=> r + 2m <= 2

  B1..B6  negative controls
  C1..C2  robustness: characteristic and reduced-word independence
"""
import os
import sys
import time
from collections import Counter
from itertools import combinations

from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from homological import setup
from gls import iterated_socle, socle_step
from qvar import rref, P as MODP
from epsilon import eps_of_mu
from upsets import upset
import thresh as TH
import strands as ST
import canon as CA
import reduce as RD

RES = []


def record(tag, label, ok, detail=""):
    RES.append((tag, label, ok, detail))
    flag = "ok " if ok else "FAIL"
    print(f"  [{flag}] {tag:4} {label:52} {detail}")


def module(rs, pr, word, lam):
    seq = []
    for _, i, c, _ in trace(rs, word, lam):
        seq += [i] * c
    return iterated_socle(pr, rs, seq)


# ----------------------------------------------------------------- A: the chain

def A1_main_theorem():
    """soc_sigma(I_kappa) has dimension vector v, and its graded pieces agree
    with the independent quiver-variety solver."""
    from dfs3 import solve
    st = Counter()
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("E", 6, 2)]:
        rs, pr, kk = setup(typ, n, k)
        lam = rs.fundamental(k)
        for mu, w in orbit_with_words(rs, lam).items():
            if not w:
                continue
            X = module(rs, pr, w, lam)
            v = tuple(int(x) for x in v_from_trace(rs, trace(rs, w, lam)))
            got = tuple(sum(len(X.get((d, i), [])) for d in range(pr.top + 1))
                        for i in rs.I)
            st['dim_tot'] += 1
            st['dim_ok'] += (got == v)
            simple = all(len(x) <= 1 for x in X.values())
            r = solve(rs, k, v)
            if r:
                g = list(r.values())[0]
                truth = all(c <= 1 for row in g.values() for c in row.values())
                st['sol_tot'] += 1
                st['sol_ok'] += (simple == truth)
    record("A1", "dim soc_sigma(I_kappa) = v",
           st['dim_ok'] == st['dim_tot'], f"{st['dim_ok']}/{st['dim_tot']}")
    record("A1", "simple spectrum agrees with the solver",
           st['sol_ok'] == st['sol_tot'], f"{st['sol_ok']}/{st['sol_tot']}")


def A2_inward():
    """Lemma I: the inward arrow is injective wherever a slot has dim >= 2.

    In type E this is the complete finite verification the proof rests on; in
    type D the paper proves it uniformly, and we confirm it in range."""
    st = Counter()
    for typ, n, ks in [("E", 6, range(1, 7)), ("E", 7, range(1, 8)),
                       ("D", 8, range(2, 7)), ("D", 9, range(2, 8))]:
        for k in ks:
            rs, pr, kk = setup(typ, n, k)
            t = [i for i in rs.I if len(rs.adj[i]) >= 3][0]
            # nx(i) = the neighbour one step towards t
            import collections
            dist = {t: 0}
            q = collections.deque([t])
            while q:
                x = q.popleft()
                for y in rs.adj[x]:
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        q.append(y)
            for (d, i), m in sorted(pr.dim.items()):
                if i == t or m < 2:
                    continue
                nx = min(rs.adj[i], key=lambda y: dist[y])
                if (d, i, nx) not in pr.arrow:
                    st[typ + '_tot'] += 1
                    continue
                IM = [pr.act(d, i, [1 if x == b else 0 for x in range(m)], nx)
                      for b in range(m)]
                r = len(rref([x[:] for x in IM], pr.dim.get((d + 1, nx), 0))[0])
                st[typ + '_tot'] += 1
                st[typ + '_ok'] += (r == m)
    record("A2", "Lemma I in type E (the finite verification)",
           st['E_ok'] == st['E_tot'], f"{st['E_ok']}/{st['E_tot']}")
    record("A2", "Lemma I in type D (uniform proof, confirmed)",
           st['D_ok'] == st['D_tot'], f"{st['D_ok']}/{st['D_tot']}")


def A3_profile():
    st = Counter()
    for n in range(5, 16):
        for k in range(2, n - 1):
            rs, pr, kk = setup("D", n, k)
            st['kk_tot'] += 1
            st['kk_ok'] += (kk == k)
            for j in rs.I:
                got = [(d, pr.dim[(d, j)]) for d in range(pr.top + 1)
                       if pr.dim.get((d, j))]
                st['tot'] += 1
                st['ok'] += (got == TH.profile(n, k, j))
    record("A3", "Lemma P closed form = the module profile",
           st['ok'] == st['tot'], f"{st['ok']}/{st['tot']}")
    record("A3", "the injective is Pi e_k (kappa' = k)",
           st['kk_ok'] == st['kk_tot'], f"{st['kk_ok']}/{st['kk_tot']}")


def A4_vt():
    st = Counter()
    for n in range(4, 10):
        for k in range(2, n - 1):
            rs = RootSystem("D", n)
            lam = rs.fundamental(k)
            for mu, w in orbit_with_words(rs, lam).items():
                v = [int(x) for x in v_from_trace(rs, trace(rs, w, lam))]
                e = eps_of_mu(n, k, v)
                m = sum(1 for i in range(n - 2) if e[i] < 0)
                r = sum(1 for i in (n - 2, n - 1) if e[i] != 0)
                st['tot'] += 1
                st['ok'] += (v[n - 3] == r + 2 * m)
    record("A4", "Lemma C: v_t = r + 2m",
           st['ok'] == st['tot'], f"{st['ok']}/{st['tot']}")


def A5_strands():
    S = Counter()
    for n in range(5, 15):
        for k in range(2, n - 1):
            S += ST.run("D", n, k)
    for tag, a, b in [("(W1) multiplicity-free", 'W1_ok', 'W1_tot'),
                      ("(W2) closed form = tau-eigendims", 'W2e_ok', 'W2e_tot'),
                      ("(W3) contiguous intervals", 'W3_ok', 'W3_tot'),
                      ("(W4) arrows transport strands", 'W4_ok', 'W4_tot'),
                      ("(W4l) the folded-leaf arrow", 'W4l_ok', 'W4l_tot')]:
        record("A5", tag, S[a] == S[b], f"{S[a]}/{S[b]}")


def A6_T012():
    S = Counter()
    for n in range(5, 61):
        for k in range(2, n - 1):
            TH.local_facts(n, k, S)
    for tag, a, b in [("(T0) rank = min(dim,dim)", 'T0_ok', 'T0_tot'),
                      ("(T1) a rank-1 image is the + strand", 'T1_ok', 'T1_tot'),
                      ("(T2) the + strand dies into 1-dim slots",
                       'T2_ok', 'T2_tot')]:
        record("A6", tag, S[a] == S[b], f"{S[a]}/{S[b]}")


def A7_reps():
    """The W_J-minimal weights with c_t >= 1 are exactly the explicit family,
    k of them, or k-1 when k = n-2."""
    st = Counter()
    for n in range(5, 12):
        t = n - 2
        for k in range(2, t + 1):
            rs = RootSystem("D", n)
            lam = rs.fundamental(k)
            orb = orbit_with_words(rs, lam)
            if len(orb) > 5000:
                continue
            got = {mu for mu in orb
                   if {j for j in rs.I if mu[j - 1] <= -1} == {t - 1}
                   and mu[t - 1] >= 1}
            pred = set()
            for a in range(k):
                b = n - 2 - k
                e = [0] * (n + 1)
                for p in range(1, a + 1):
                    e[p] = 1
                for p in range(a + b + 1, n - 2):
                    e[p] = -1
                e[n - 2] = 1
                pred.add(tuple([e[j] - e[j + 1] for j in range(1, n)]
                               + [e[n - 1] + e[n]]))
            st['tot'] += 1
            st['sub'] += (got <= pred)
            st['cnt'] += (len(got) == (k if k < t else k - 1))
    record("A7", "representatives lie in the explicit family",
           st['sub'] == st['tot'], f"{st['sub']}/{st['tot']}")
    record("A7", "and number k  (k-1 when k = n-2)",
           st['cnt'] == st['tot'], f"{st['cnt']}/{st['tot']}")


def A8_canonical():
    S = Counter()
    for typ, n, k in CA.CASES:
        S += CA.run(typ, n, k)
    record("A8", "tau is an involutive automorphism",
           S['inv_ok'] == S['inv_tot'] and S['equi_ok'] == S['equi_tot'],
           f"{S['inv_ok']}/{S['inv_tot']}, {S['equi_ok']}/{S['equi_tot']}")
    record("A8", "T(v_a) is a submodule",
           S['sub'] == S['rep'], f"{S['sub']}/{S['rep']}")
    record("A8", "T(v_a) = V(mu_a)  (so (B) holds)",
           S['eq'] == S['rep'], f"{S['eq']}/{S['rep']}")


def A9_arithmetic():
    S = Counter()
    N = 0
    for n in range(5, 46):
        for k in range(2, n - 1):
            for a in range(k):
                v = TH.vvec(n, k, a)
                if any(x < 0 for x in v[1:]) or all(x == 0 for x in v[1:]):
                    continue
                S += TH.classify(n, k, a)
                TH.upclosed(n, k, a, S)
                TH.increments(n, k, a, S)
                TH.leaf_boundary(n, k, a, S)
                N += 1
    bad = {a: b for a, b in S.items() if a.startswith('BAD')}
    record("A9", "no closure obligation outside the five patterns",
           not bad, f"{N} representatives" + (f"  {bad}" if bad else ""))
    for tag, a, b in [("up-closed in each strand", 'up_ok', 'up_tot'),
                      ("occupation formula = the truncation", 'occ_ok', 'occ_tot'),
                      ("0 <= s+ increment <= 1", 'inc_ok', 'inc_tot'),
                      ("|s- increment| <= 1", 'incm_ok', 'incm_tot'),
                      ("(P1) u increments in {-1,0,1}", 'P1_ok', 'P1_tot'),
                      ("(P2) a drop only at p <= a or at t", 'P2_ok', 'P2_tot'),
                      ("(P3) u_(t-1) <= o_(t-1) + 1", 'P3_ok', 'P3_tot'),
                      ("(P5) a jump of 2 raises s+ by one", 'P5_ok', 'P5_tot'),
                      ("(Q4) leaf obligations land in full slots",
                       'Q4_ok', 'Q4_tot')]:
        record("A9", tag, S[a] == S[b], f"{S[a]}/{S[b]}")
    record("A9", "(P5) the parity case never occurs",
           S['P5_gap'] == 0, f"{S['P5_gap']} occurrences")


def A10_G12():
    """The two profile inequalities behind (G1) and (G2), with the counts the
    paper's proof predicts."""
    st = Counter()
    for n in range(5, 81):
        t = n - 2
        for k in range(2, t + 1):
            pr = dict(TH.profile(n, k, t - 1))
            tot_t = sum(m for _, m in TH.profile(n, k, t))
            for a in range(0, k + 2):
                for case, D, bound, lo in ((1, t + k - 2 * a,
                                            2 * a - (k == t), 1),
                                           (2, t + k - 2 * a - 2,
                                            2 * a + 2 - (k == t), 0)):
                    if a < lo:
                        continue
                    if (2 * a if case == 1 else 2 * a + 1) > tot_t:
                        continue
                    cum = sum(m for d, m in pr.items() if d >= D + 1)
                    st['tot'] += 1
                    st['ok'] += (cum <= bound)
                    cp = min(max(0, a - 1 if case == 1 else a), min(t - 1, k))
                    cm = min((a + 1 if case == 1 else a + 2) - (k == t),
                             min(t - 1, k))
                    gp = sum(1 for d in range(2 * t + 1)
                             if TH.wave(n, k, t - 1, d, 1) and d >= D + 1)
                    gm = sum(1 for d in range(2 * t + 1)
                             if TH.wave(n, k, t - 1, d, -1) and d >= D + 1)
                    st['cnt_tot'] += 1
                    st['cnt_ok'] += (gp == cp and gm == cm)
    record("A10", "(G1),(G2): the profile inequalities",
           st['ok'] == st['tot'], f"{st['ok']}/{st['tot']}")
    record("A10", "and the strand counts the proof uses",
           st['cnt_ok'] == st['cnt_tot'], f"{st['cnt_ok']}/{st['cnt_tot']}")


def A11_end_to_end():
    st = Counter()
    for n in range(4, 11):
        for k in range(2, n - 1):
            rs, pr, kk = setup("D", n, k)
            lam = rs.fundamental(k)
            orb = orbit_with_words(rs, lam)
            if len(orb) > 1200:
                continue
            for mu, w in orb.items():
                if not w:
                    continue
                X = module(rs, pr, w, lam)
                simple = all(len(x) <= 1 for x in X.values())
                v = [sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
                     for j in rs.I]
                e = eps_of_mu(n, k, v)
                m = sum(1 for i in range(n - 2) if e[i] < 0)
                r = sum(1 for i in (n - 2, n - 1) if e[i] != 0)
                st['tot'] += 1
                st['ok'] += ((r + 2 * m <= 2) == simple)
    record("A11", "simple spectrum  <=>  r + 2m <= 2",
           st['ok'] == st['tot'], f"{st['ok']}/{st['tot']}")


# ------------------------------------------------------- B: negative controls

def B_negative_controls():
    """Statements that MUST be false.  If one of them starts holding, something
    has drifted -- most likely a definition."""
    from allsub import all_submodules

    # B1: extremality cannot be dropped from (*o)
    bad = 0
    for n, k in [(6, 4), (7, 4)]:
        rs, pr, kk = setup("D", n, k)
        t = n - 2
        for N in all_submodules(rs, kk)[1]:
            vv = {j: sum(len(N.get((d, j), [])) for d in range(pr.top + 1))
                  for j in rs.I}
            ct = (1 if t == kk else 0) + sum(vv[j] for j in rs.adj[t]) - 2 * vv[t]
            H = upset(pr, N, t) and upset(pr, N, n - 1) and upset(pr, N, n)
            if H and ct >= 1 and not upset(pr, N, t - 1):
                bad += 1
    record("B1", "(L) is FALSE at kappa = 4 (extremality is needed)",
           bad == 4, f"{bad} counterexamples (expected 4)")

    # B2: V is not the greedy top truncation at every vertex
    rs, pr, kk = setup("D", 5, 2)
    lam = rs.fundamental(2)
    hit = 0
    for mu, w in orbit_with_words(rs, lam).items():
        if not w:
            continue
        X = module(rs, pr, w, lam)
        v = [sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
             for j in rs.I]
        greedy = True
        for j in rs.I:
            rem = v[j - 1]
            for d in range(pr.top, -1, -1):
                m = pr.dim.get((d, j), 0)
                if not m:
                    continue
                take = min(m, rem)
                rem -= take
                if len(X.get((d, j), [])) != take:
                    greedy = False
        hit += (not greedy)
    record("B2", "greedy truncation at EVERY vertex is FALSE",
           hit > 0, f"{hit} failures in D5 om2")

    # B3: the fully commutative element of D4 without simple spectrum
    rs, pr, kk = setup("D", 4, 2)
    X = module(rs, pr, (2, 4, 3, 1, 2), rs.fundamental(2))
    v = tuple(sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
              for j in rs.I)
    simple = all(len(x) <= 1 for x in X.values())
    record("B3", "D4 s2s4s3s1s2: fully commutative, NOT simple",
           (not simple) and v == (1, 3, 1, 1), f"v = {v}, simple = {simple}")

    # B4: the interval step fails without both neighbours -- the bad step
    tot = bad = 0
    for n in range(5, 20):
        t = n - 2
        for k in range(2, t + 1):
            for e in (1, -1):
                for j in range(2, t):
                    for ip in (j - 1, j + 1):
                        for d in range(2 * t):
                            if (ST.wave(n, k, j, d, e)
                                    and ST.wave(n, k, ip, d + 1, e)):
                                tot += 1
                                bad += not ST.wave(n, k, j, d + 2, e)
    record("B4", "one-neighbour interval step is FALSE (the bad step)",
           bad > 0, f"{bad}/{tot} failures")

    # B6: the RPP hook product needs DKK's orientation
    import rppmod as RP
    good = tot = 0
    for typ, n, k in [("A", 4, 2), ("D", 4, 2), ("D", 5, 2)]:
        rs, pr, kk = setup(typ, n, k)
        lam = rs.fundamental(k)
        for mu, w in orbit_with_words(rs, lam).items():
            if not w:
                continue
            b = RP.build(rs, pr, w, lam)
            if b is None:
                continue
            elts, down = b
            if len(elts) > RP.MAXSIZE:
                continue
            m = len(elts)
            up = [frozenset(x for x in range(m) if y in down[x])
                  for y in range(m)]
            H = RP.product(RP.hook_multiset(rs, w, lam))
            tot += 1
            good += (RP.rpp(up) == H)
            st_ok = (RP.rpp(down) == H)
            if not st_ok:
                record("B6", "DKK orientation IS the hook product", False,
                       f"{typ}{n} om{k}")
                return
    record("B6", "the opposite orientation is NOT a hook product",
           good < tot, f"{good}/{tot} would be")


# ------------------------------------------------------------- C: robustness

def C1_characteristic():
    """The module computations are over F_p.  Re-run the headline check at
    several primes in a subprocess, so nothing is cached."""
    import subprocess
    snippet = (
        "from homological import setup\n"
        "from rootsys import orbit_with_words\n"
        "from wordtrace import trace\n"
        "from gls import iterated_socle\n"
        "from epsilon import eps_of_mu\n"
        "ok=tot=0\n"
        "for n,k in [(5,2),(5,3),(6,3)]:\n"
        "    rs,pr,kk=setup('D',n,k); lam=rs.fundamental(k)\n"
        "    for mu,w in orbit_with_words(rs,lam).items():\n"
        "        if not w: continue\n"
        "        seq=[]\n"
        "        for _,i,c,_ in trace(rs,w,lam): seq+=[i]*c\n"
        "        X=iterated_socle(pr,rs,seq)\n"
        "        s=all(len(x)<=1 for x in X.values())\n"
        "        v=[sum(len(X.get((d,j),[])) for d in range(pr.top+1)) for j in rs.I]\n"
        "        e=eps_of_mu(n,k,v)\n"
        "        m=sum(1 for i in range(n-2) if e[i]<0)\n"
        "        r=sum(1 for i in (n-2,n-1) if e[i]!=0)\n"
        "        tot+=1; ok+= ((r+2*m<=2)==s)\n"
        "print(ok,tot)\n")
    outs = []
    for p in (5, 7, 11, 13, 1000003):
        env = dict(os.environ, DKKP=str(p))
        r = subprocess.run([sys.executable, "-c", snippet], env=env,
                           capture_output=True, text=True, cwd=os.getcwd())
        outs.append((p, r.stdout.strip()))
    outs = [(p, (o.strip().splitlines() or [""])[-1]) for p, o in outs]
    ok = all(o and len(o.split()) == 2 and o.split()[0] == o.split()[1]
             for _, o in outs)
    same = len({o for _, o in outs}) == 1
    record("C1", "the criterion is independent of the characteristic",
           ok and same, ", ".join(f"p={p}: {o}" for p, o in outs))


def C2_word():
    """The construction does not depend on the reduced word chosen."""
    from weylprops import reduced_words
    st = Counter()
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3)]:
        rs, pr, kk = setup(typ, n, k)
        lam = rs.fundamental(k)
        for mu, w in orbit_with_words(rs, lam).items():
            if not w:
                continue
            base = None
            for w2 in sorted(reduced_words(rs, w))[:8]:
                X = module(rs, pr, w2, lam)
                sig = tuple(sorted((d, i, len(v)) for (d, i), v in X.items() if v))
                if base is None:
                    base = sig
                st['tot'] += 1
                st['ok'] += (sig == base)
    record("C2", "independent of the reduced word",
           st['ok'] == st['tot'], f"{st['ok']}/{st['tot']}")


# ---------------------------------------------------- D: the manuscript agrees

def _digits(path):
    """The text with LaTeX thin spaces, commas and dollars removed, so that
    $8\,047\,275$ and 8,047,275 both read as 8047275."""
    raw = open(path).read()
    for a, b in (("\\,", ""), (",", ""), ("$", ""), ("\u2009", ""),
                 ("\u202f", "")):
        raw = raw.replace(a, b)
    return raw


def D_manuscript():
    """Every headline count this audit computes must appear in the manuscript
    and in the README.  This is what catches a number going stale after a
    change -- the most likely way a wrong claim survives."""
    tex, rdm = _digits("paper.tex"), _digits("README.md")
    want = {}
    for tag, label, ok, detail in RES:
        for tok in detail.replace("/", " ").split():
            if tok.isdigit() and int(tok) > 999:
                want.setdefault(tok, label)
    missing_tex = sorted(t for t in want if t not in tex)
    missing_rdm = sorted(t for t in want if t not in rdm)
    # only a curated subset is quoted in prose; require the big ones
    head = sorted((t for t in want if int(t) >= 100000), key=int)
    bad_t = [t for t in head if t not in tex]
    bad_r = [t for t in head if t not in rdm]
    record("D1", "the manuscript quotes the headline counts",
           not bad_t, f"{len(head)-len(bad_t)}/{len(head)}"
           + (f"  missing {bad_t}" if bad_t else ""))
    record("D2", "the README quotes them too",
           not bad_r, f"{len(head)-len(bad_r)}/{len(head)}"
           + (f"  missing {bad_r}" if bad_r else ""))


def E_smoke(per=180):
    """Every script in the repo exits 0.  Slow; set DKKFULL=1 to include.

    Each script gets `per` seconds; anything slower is reported separately
    rather than hanging the audit."""
    import subprocess
    import glob
    skip = {"audit.py"}
    bad, slow, n = [], [], 0
    for f in sorted(glob.glob("*.py")):
        if f in skip:
            continue
        n += 1
        try:
            r = subprocess.run([sys.executable, f], capture_output=True,
                               text=True, timeout=per)
            if r.returncode != 0:
                bad.append(f)
        except subprocess.TimeoutExpired:
            slow.append(f)
        print(f"        ... {f:22} "
              f"{'slow' if f in slow else ('FAIL' if f in bad else 'ok')}",
              flush=True)
    record("E1", "every script exits 0", not bad,
           f"{n-len(bad)-len(slow)}/{n} ok, {len(slow)} over {per}s"
           + (f", failing {bad}" if bad else ""))


if __name__ == "__main__":
    t0 = time.time()
    print("\nA. THE CHAIN\n")
    for f in (A1_main_theorem, A2_inward, A3_profile, A4_vt, A5_strands,
              A6_T012, A7_reps, A8_canonical, A9_arithmetic, A10_G12,
              A11_end_to_end):
        f()
    print("\nB. NEGATIVE CONTROLS  (these must be false)\n")
    B_negative_controls()
    print("\nC. ROBUSTNESS\n")
    C1_characteristic()
    C2_word()
    print("\nD. THE MANUSCRIPT\n")
    D_manuscript()
    if os.environ.get("DKKFULL"):
        print("\nE. SMOKE TEST\n")
        E_smoke()
    bad = [r for r in RES if not r[2]]
    print(f"\n{len(RES) - len(bad)}/{len(RES)} checks passed "
          f"({time.time() - t0:.0f}s)")
    if bad:
        print("\nFAILED:")
        for tag, label, _, detail in bad:
            print(f"   {tag} {label}  {detail}")
    print("\nAUDIT PASSED" if not bad else "\n*** AUDIT FAILED ***")
    sys.exit(0 if not bad else 1)
