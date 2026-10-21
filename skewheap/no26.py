"""Naruse-Okada Definition 2.6 (recursive hook monomial), univariate, implemented from the
quoted definition and validated against the hooks extracted from the RPP series.

Def 2.6 (as quoted): if u is not the top of any d_k-interval, z[H_P(u)] = prod_{w <= u}
z_{c(w)}; if u is the top of a d_k-interval [v,u] with sides x and y, then
z[H_P(u)] = z[H_P(x)] z[H_P(y)] / z[H_P(v)].
Setting every z_i = q turns products into sums of exponents, so

    h(u) = #{w : w <= u}                     if u tops no d_k-interval,
    h(u) = h(x) + h(y) - h(v)                if [v,u] is a d_k-interval with sides x,y.

Orientation here is the LITERATURE's: v < u means v is below u, and a d_k-interval [v,u]
has v as its unique minimum and u as its unique maximum.
usage: python3 no26.py [MAXLEN=9] [PMAX=14] [TYPES=A4,A5,D4,D5]
"""
import sys, collections
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from crit import is_product
from novel import dual, relabel_natural, ser

def leq_matrix(down):
    """leq[a][b] True iff a <= b.  down[u] = elements immediately below u."""
    n = len(down)
    leq = [[False]*n for _ in range(n)]
    for a in range(n): leq[a][a] = True
    for _ in range(n):
        for u in range(n):
            for v in down[u]:
                for a in range(n):
                    if leq[a][v] and not leq[a][u]: leq[a][u] = True
    return leq

def dk_intervals(down, leq):
    """all (v,u,x,y) with [v,u] a d_k-interval: unique min v, unique max u, exactly two
    incomparable elements x,y, and the rest a chain below both plus a chain above both."""
    n = len(down); out = []
    for v in range(n):
        for u in range(n):
            if v == u or not leq[v][u]: continue
            I = [w for w in range(n) if leq[v][w] and leq[w][u]]
            if len(I) < 4: continue
            inc = [(a, b) for i, a in enumerate(I) for b in I[i+1:]
                   if not leq[a][b] and not leq[b][a]]
            if len(inc) != 1: continue
            x, y = inc[0]
            L = [w for w in I if w not in (x, y) and leq[w][x] and leq[w][y]]
            U = [w for w in I if w not in (x, y) and leq[x][w] and leq[y][w]]
            if len(L) + len(U) + 2 != len(I): continue
            # d_k(1) is the double-tailed diamond with 2k-2 elements: a chain of k-2
            # below the incomparable pair and a chain of k-2 above it.  So |L| = |U|.
            if len(L) != len(U): continue
            chain = lambda X: all(leq[a][b] or leq[b][a] for a in X for b in X)
            if not (chain(L) and chain(U)): continue
            if v not in L or u not in U: continue
            out.append((v, u, x, y))
    return out

def no_hooks(down):
    """h(u) for all u, by Def 2.6; returns None if the recursion is ill-posed."""
    n = len(down); leq = leq_matrix(down)
    tops = {}
    for (v, u, x, y) in dk_intervals(down, leq):
        tops.setdefault(u, []).append((v, x, y))
    order = sorted(range(n), key=lambda u: sum(1 for w in range(n) if leq[w][u]))
    h = {}
    for u in order:
        if u in tops:
            v, x, y = tops[u][0]
            if x not in h or y not in h or v not in h: return None
            h[u] = h[x] + h[y] - h[v]
        else:
            h[u] = sum(1 for w in range(n) if leq[w][u])
    return [h[u] for u in range(n)]

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "A4,A5,D4,D5").split(",")
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
            if n < 2 or n > PMAX: continue
            dn = relabel_natural(down)
            dc = (len(dcomplete(relabel_natural(dual(down)))) == 0)
            try: ok, b, m = is_product(ser(dn), 29, n)
            except AssertionError: continue
            prod = bool(ok) and sum(m.values()) == n
            if not prod: continue
            series_hooks = sorted(k for k, v in m.items() for _ in range(v))
            # Def 2.6 lives in the convention where the "corner" is the unique MAXIMUM,
            # i.e. the dual of this codebase's heap orientation -- the same dual that
            # dcomplete() requires.  (Checked on the Young diagram (2,1): hooks 3,1,1.)
            nh = no_hooks(relabel_natural(dual(dn)))
            if nh is None: stat[(dc, 'recursion ill-posed')] += 1; continue
            agree = sorted(nh) == series_hooks
            stat[(dc, 'Def 2.6 = series hooks' if agree else 'DISAGREE')] += 1
            if not agree and len(ex) < 4:
                ex.append((t, w, n, dc, series_hooks, sorted(nh)))
    print(f"types {TY}, words to length {ML}, |P| <= {PMAX}")
    print("product heaps only; comparing Def 2.6's recursive hooks with the hooks read off")
    print("the RPP series\n")
    for k in sorted(stat, key=str):
        print(f"   d-complete={str(k[0]):5} {k[1]:24} {stat[k]}")
    for a in ex:
        print(f"\n   DISAGREE {a[0]} w={a[1]} |P|={a[2]} d-complete={a[3]}")
        print(f"      series hooks  = {a[4]}")
        print(f"      Def 2.6 hooks = {a[5]}")

if __name__ == "__main__":
    main()
