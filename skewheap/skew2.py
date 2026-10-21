"""Naruse-Okada Theorem 1.2, with excited PEAKS.

    sum_{pi in A(P\\F)} q^{|pi|}  =  sum_D  ( prod_{v in B(D)} q^{h(v)} )
                                            / ( prod_{v in P\\D} (1 - q^{h(v)}) )

D ranges over excited diagrams of the order filter F; B(D) is the set of excited peaks.
Excitation (Def 3.1 + the peak rule): u in D is D-active if some d_k-interval [v,u] has
v not in D and every z in [v,u] that is covered by u or covers v lies outside D.  Then
    alpha_u(D) = (D \\ {u}) union {v},
    B(alpha_u(D)) = ( B(D) \\ { z in [v,u] : z covered by u or covers v } ) union {PEAK},
where PEAK is the element the transcription leaves ambiguous.  We try both readings and
let the theorem itself decide, on d-complete posets where it is proved.
usage: python3 skew2.py [MAXLEN=7] [PMAX=9] [TYPES=A4] [MAXF=2] [PEAK=u|v]
"""
import sys, itertools, collections
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from novel import dual, relabel_natural, ser
from no26 import leq_matrix, dk_intervals, no_hooks
from filterhook import hooks_of
from skew import to_proctor, series_prod

def heap_sub_series(down, heap_ids, N):
    """(induced subposet)-partition series, computed in the HEAP orientation that ser()
    was validated against -- feeding ser() a Proctor-orientation poset gives the series of
    the order-dual, which is a different generating function."""
    S = sorted(heap_ids); pos = {u: i for i, u in enumerate(S)}
    sp = [frozenset(pos[v] for v in down[u] if v in pos) for u in S]
    return ser(relabel_natural(sp), N)

def build(dP):
    n = len(dP); leq = leq_matrix(dP)
    cov = [[v for v in range(n) if v != u and leq[v][u] and
            not any(w not in (u, v) and leq[v][w] and leq[w][u] for w in range(n))]
           for u in range(n)]           # cov[u] = elements covered by u
    dks = {}
    for (v, u, x, y) in dk_intervals(dP, leq):
        I = [w for w in range(n) if leq[v][w] and leq[w][u]]
        # the "rim" of the interval: covered by u, or covering v
        rim = set(w for w in I if w in cov[u]) | set(w for w in I if u != w and v in cov[w])
        dks.setdefault(u, []).append((v, I, frozenset(rim)))
    return leq, dks

def excited(dP, F, peak, cap=40000):
    leq, dks = build(dP)
    start = (frozenset(F), frozenset())
    seen = {start}; stack = [start]
    while stack:
        D, B = stack.pop()
        for u in list(D):
            for (v, I, rim) in dks.get(u, []):
                if v in D: continue
                if any(z in D for z in rim): continue
                D2 = frozenset((set(D) - {u}) | {v})
                B2 = frozenset((set(B) - set(rim)) | {u if peak == 'u' else v})
                s = (D2, B2)
                if s not in seen:
                    seen.add(s); stack.append(s)
                    if len(seen) > cap: return None
    return seen

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "A4").split(",")
    MAXF = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    PEAK = sys.argv[5] if len(sys.argv) > 5 else 'u'
    N = 18
    stat = collections.Counter(); ex = []
    for t in TY:
        typ, rk = t[0], int(t[1:]); C = cartan(typ, rk)
        out = []
        def rec(w):
            if w and is_reduced(C, w): out.append(tuple(w))
            if len(w) == ML: return
            for i in range(rk): rec(w + [i])
        rec([])
        for w in out:
            down = heap_down(C, list(w)); n = len(down)
            if n < 3 or n > PMAX: continue
            dc = (len(dcomplete(relabel_natural(dual(down)))) == 0)
            hp = (hooks_of(relabel_natural(down)) is not None)   # hook length poset?
            dP, col, order = to_proctor(down, list(w))
            h = no_hooks(dP)
            if h is None: continue
            leq = leq_matrix(dP)
            ups = [frozenset(x for x in range(n) if leq[u][x]) for u in range(n)]
            filters = set()
            for r in range(1, MAXF+1):
                for S in itertools.combinations(range(n), r):
                    Fl = frozenset().union(*[ups[u] for u in S])
                    if len(Fl) <= MAXF+2: filters.add(Fl)
            for Fl in filters:
                if not Fl or len(Fl) >= n: continue
                E = excited(dP, Fl, PEAK)
                if E is None: stat['too many states'] += 1; continue
                # is B a function of D?
                byD = collections.defaultdict(set)
                for D, B in E: byD[D].add(B)
                if any(len(v) > 1 for v in byD.values()): stat['B not a function of D'] += 1
                lhs = heap_sub_series(down, [order[u] for u in range(n) if u not in Fl], N)
                rhs = [0]*(N+1)
                for D, B in E:
                    shift = sum(h[v] for v in B)
                    if shift > N: continue
                    pr = series_prod([h[u] for u in range(n) if u not in D], N)
                    for d in range(N+1-shift): rhs[d+shift] += pr[d]
                agree = lhs[:N+1] == rhs[:N+1]
                stat[('d-complete' if dc else ('hook-length, not d-complete' if hp
                      else 'not a hook length poset'), agree)] += 1
                if not agree and hp and not dc and len(ex) < 4:
                    ex.append((t, w, n, sorted(Fl), len(E), lhs[:7], rhs[:7]))
    print(f"types {TY}, len<={ML}, |P|<={PMAX}, filters<={MAXF+2}, peak rule adds '{PEAK}'\n")
    for k in sorted(stat, key=str): print(f"   {k} -> {stat[k]}")
    for a in ex:
        print(f"\n   HOOK-LENGTH (not d-complete) MISMATCH {a[0]} w={a[1]} |P|={a[2]} "
              f"F={a[3]} #states={a[4]}")
        print(f"      LHS = {a[5]}")
        print(f"      RHS = {a[6]}")

if __name__ == "__main__":
    main()
