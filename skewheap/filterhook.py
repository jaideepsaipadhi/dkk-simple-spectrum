"""Per-element hooks from series data alone.

For a Young diagram, the hook of a cell u equals the hook of the CORNER of the subdiagram
sitting weakly above u, and in a product poset the corner carries the LARGEST hook.  So
the candidate rule is
        h_P(u)  =  max( hook multiset of the principal order filter { w >= u } ),
which needs only the RPP series of each principal filter -- no d-completeness, no Def 2.6.
Validate against Def 2.6 on d-complete heaps, then use it as ground truth elsewhere.
usage: python3 filterhook.py [MAXLEN=8] [PMAX=12] [TYPES=A4,A5,D4]
"""
import sys, collections
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from crit import is_product
from novel import dual, relabel_natural, ser
from no26 import leq_matrix, no_hooks

def hooks_of(d):
    n = len(d)
    try: ok, b, m = is_product(ser(d), 29, n)
    except AssertionError: return None
    if not (ok and sum(m.values()) == n): return None
    return sorted(k for k, v in m.items() for _ in range(v))

def filter_hooks(d):
    """h(u) = max hook of the principal filter {w >= u}, for every u; None if some filter
    is not a product poset."""
    n = len(d); leq = leq_matrix(d)
    out = []
    for u in range(n):
        S = sorted(w for w in range(n) if leq[u][w])     # w >= u
        pos = {w: i for i, w in enumerate(S)}
        subp = [frozenset(pos[v] for v in d[w] if v in pos) for w in S]
        h = hooks_of(relabel_natural(subp))
        if h is None: return None
        out.append(max(h))
    return out

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "A4,A5,D4").split(",")
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
            th = hooks_of(dn)
            if th is None: continue
            fh = filter_hooks(dn)
            if fh is None: stat[(dc, 'a filter is not a product')] += 1; continue
            agree_multiset = sorted(fh) == th
            nh = no_hooks(relabel_natural(dual(dn)))
            agree_no26 = (nh is not None and sorted(nh) == sorted(fh))
            stat[(dc, f'filter-rule matches series: {agree_multiset}')] += 1
            if dc: stat[(dc, f'filter-rule = Def 2.6: {agree_no26}')] += 1
            if not agree_multiset and len(ex) < 4:
                ex.append((t, w, n, dc, th, sorted(fh)))
    print(f"types {TY}, words to length {ML}, |P| <= {PMAX}")
    print("product heaps; candidate rule h(u) = max hook of the filter above u\n")
    for k in sorted(stat, key=str):
        print(f"   d-complete={str(k[0]):5} {k[1]:34} {stat[k]}")
    for a in ex:
        print(f"\n   MISMATCH {a[0]} w={a[1]} |P|={a[2]} d-complete={a[3]}")
        print(f"      series hooks      = {a[4]}")
        print(f"      filter-rule hooks = {a[5]}")

if __name__ == "__main__":
    main()
