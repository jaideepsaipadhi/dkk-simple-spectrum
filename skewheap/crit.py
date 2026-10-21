"""CRITERION P: a closed-form test for "RPP_P is a product", from Lemma U.

Lemma U's identity N*[q^N] log RPP_P = sum_{d|N} d*m_d inverts by Mobius:

    d * m_d = sum_{e | d} mu(d/e) * a_e,        a_e := e * [q^e] log RPP_P.

The a_e are determined by the series.  So the candidate multiplicities m_d are EXPLICIT, and

    RPP_P is a product of 1/(1-q^{h_i})  <=>  every m_d is a nonnegative integer

(the count sum_d m_d is then automatically |P|, by pole order at q=1).  No greedy division, no
candidate enumeration: a formula, with the failure located at a specific d.

usage: python3 crit.py            (agreement check against the greedy test, all posets |P|<=5)
       python3 crit.py cases      (the m_d of the two-unclosed-diamond classes)"""
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from sympy import mobius, divisors

def logseries(f, N):
    """[q^0..q^{N-1}] of log(f), f a coefficient list with f[0] == 1."""
    f = [Fr(v) for v in f[:N]]
    assert f[0] == 1
    L = [Fr(0)]*N
    # f' / f  then integrate:  (log f)' = f'/f
    # do it directly: L[n] = f[n] - (1/n) * sum_{k=1}^{n-1} k*L[k]*f[n-k]
    for n in range(1, N):
        s = Fr(0)
        for k in range(1, n):
            s += k*L[k]*f[n-k]
        L[n] = f[n] - s/n
    return L

def mults(f, N):
    """candidate multiplicities m_d for d < N; returns dict d -> Fraction."""
    L = logseries(f, N)
    a = {e: e*L[e] for e in range(1, N)}
    out = {}
    for d in range(1, N):
        s = Fr(0)
        for e in divisors(d):
            s += mobius(d//e)*a[e]
        out[d] = s/d
    return out

def is_product(f, N, n):
    m = mults(f, N)
    bad = [(d, m[d]) for d in sorted(m) if m[d] < 0 or m[d].denominator != 1]
    if bad: return False, bad[0], m
    return True, None, m

def main():
    from orient import rpp
    if len(sys.argv) > 1 and sys.argv[1] == 'cases':
        from nsl import cartan, heap_down
        cases = [('A4', (0,2,1,0,3,2,1,3)), ('A4', (0,1,0,3,2,1,0,3,2)),
                 ('A4', (1,2,1,0,3,2,1,0,2)), ('A4', (1,0,2,3,2,1,0,2,1))]
        for t, w in cases:
            C = cartan(t[0], int(t[1:])); down = heap_down(C, list(w)); n = len(down)
            N = 40
            ok, bad, m = is_product(list(rpp(down, N)), N, n)
            nz = {d: m[d] for d in sorted(m) if m[d] != 0 and d <= 20}
            print(f"{t} {w} n={n}: product={ok}")
            print(f"    first bad d: {bad}")
            print(f"    m_d (d<=20, nonzero): {nz}")
        return
    from allposets import posets
    from fractions import Fraction
    def greedy(down, N, n):
        f = [Fr(v) for v in list(rpp(down, N))[:N]]; hs = []
        for d in range(1, N):
            c = f[d]
            if c < 0 or c != int(c): return False
            for _ in range(int(c)):
                hs.append(d)
                nf = f[:]
                for i in range(N-1, -1, -1):
                    if i-d >= 0: nf[i] = f[i]-f[i-d]
                f = nf
                if len(hs) > n: return False
        return all(v == 0 for v in f[1:])
    tot = dis = 0
    for n in range(1, 6):
        for down in posets(n):
            N = 36
            tot += 1
            a = greedy(down, N, n)
            b, _, _ = is_product(list(rpp(down, N)), N, n)
            if a != b:
                dis += 1
                if dis < 4: print('DISAGREE', n, down, 'greedy', a, 'criterion', b)
    print(f"agreement check: {tot} posets, {dis} disagreements")

if __name__ == "__main__":
    main()
