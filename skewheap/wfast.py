"""Fast W_P by dynamic programming over (mask, last element).

Enumerating linear extensions one by one is hopeless past |P| ~ 12.  But
comaj(w) = sum over descent positions i of (n - i), and whether position i is a descent
depends only on w[i-1] and w[i], while the weight (n - i) depends only on i = popcount(mask).
So a DP over states (set of elements already placed, last element placed) computes the whole
polynomial W_P, with 2^n * n states.

Exports: Wvec_fast(down) -> coefficient list; elinext(down) -> number of linear extensions;
degW(down) -> deg W_P.  Validated against the linear-extension enumeration in s3a9.py.
usage: python3 wfast.py        (self-check against the slow version)"""
import sys, collections
sys.path.insert(0, '.')

def Wvec_fast(down):
    n = len(down)
    need = [0]*n
    for u in range(n):
        for v in down[u]: need[u] |= 1 << v
    # dp[(mask, last)] = Counter(degree -> count)
    cur = {}
    for u in range(n):
        if need[u] == 0:
            cur[(1 << u, u)] = collections.Counter({0: 1})
    for step in range(1, n):
        nxt = {}
        for (mask, last), poly in cur.items():
            for u in range(n):
                if mask >> u & 1: continue
                if need[u] & ~mask: continue
                add = (n - step) if last > u else 0
                key = (mask | (1 << u), u)
                tgt = nxt.setdefault(key, collections.Counter())
                for d, c in poly.items(): tgt[d + add] += c
        cur = nxt
    tot = collections.Counter()
    for poly in cur.values():
        for d, c in poly.items(): tot[d] += c
    out = [0]*(max(tot)+1) if tot else [1]
    for d, c in tot.items(): out[d] = c
    return out

def elinext(down):
    return sum(Wvec_fast(down))

def degW(down):
    return len(Wvec_fast(down)) - 1

def main():
    from nsl import cartan, heap_down, is_reduced
    from s3a9 import Wvec
    import itertools
    bad = tot = 0
    for typ, rk in (('A', 4), ('D', 4), ('B', 3)):
        C = cartan(typ, rk)
        words = []
        def grow(w):
            if w: words.append(tuple(w))
            if len(w) == 7: return
            for i in range(rk):
                if is_reduced(C, w+[i]): grow(w+[i])
        grow([])
        for w in words:
            down = heap_down(C, list(w))
            if len(down) > 9: continue
            tot += 1
            if Wvec_fast(down) != Wvec(down):
                bad += 1
                if bad < 4: print('MISMATCH', typ, rk, w)
    print(f"self-check: {tot} heaps, {bad} mismatches")

if __name__ == "__main__":
    main()
