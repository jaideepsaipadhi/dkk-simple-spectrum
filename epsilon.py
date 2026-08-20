"""
Type D in epsilon-coordinates.

Bourbaki D_n: alpha_i = e_i - e_{i+1} (i <= n-1), alpha_n = e_{n-1} + e_n.
The trivalent node is t = n-2.  For lambda = omega_k = e_1 + ... + e_k
(k <= n-2), every mu in W.lambda is  mu = SUM_{i in S} eta_i e_i  with |S| = k
and signs eta_i = +-1.

Writing lambda - mu = SUM v_i alpha_i and comparing coefficients of e_j:
    v_1 = c_1,   v_j = c_1 + ... + c_j   (j <= n-2)
so
    v_t = v_{n-2} = SUM_{j=1}^{n-2} c_j
        = k - SUM_{j in S, j <= n-2} eta_j.

With r = |S cap {n-1, n}| and m = #{j in S, j <= n-2 : eta_j = -1}:

        ****  v_t = r + 2m  ****

Test this identity.
"""
import json
from rootsys import RootSystem
from wordtrace import trace, v_from_trace
from rootsys import orbit_with_words

def eps_of_mu(n, k, v):
    """recover mu in epsilon coordinates from v (= root coords of lambda-mu)"""
    c = [0]*(n+1)                       # 1-indexed
    c[1] = v[0]
    for j in range(2, n-1):
        c[j] = v[j-1] - v[j-2]
    c[n-1] = v[n-2] - v[n-3] + v[n-1]
    c[n]   = v[n-1] - v[n-2]
    lam = [0]*(n+1)
    for j in range(1, k+1): lam[j] = 1
    return [lam[j] - c[j] for j in range(1, n+1)]

bad = tot = 0
ex = []
for n in range(4, 9):
    rs = RootSystem("D", n); t = n-2
    for k in range(1, n-1):                     # non-spin fundamental weights
        lam = rs.fundamental(k)
        for mu, word in orbit_with_words(rs, lam).items():
            if not word: continue
            v = [int(x) for x in v_from_trace(rs, trace(rs, word, lam))]
            e = eps_of_mu(n, k, v)
            if any(abs(x) not in (0,1) for x in e) or sum(abs(x) for x in e) != k:
                bad += 1
                if len(ex)<3: ex.append(("bad eps", n, k, v, e)); continue
            r = sum(1 for i in (n-2, n-1) if e[i] != 0)     # 0-indexed: e[n-2]=e_{n-1}
            m = sum(1 for i in range(0, n-2) if e[i] == -1)
            tot += 1
            if v[t-1] != r + 2*m:
                bad += 1
                if len(ex)<6: ex.append((n,k,v,e,r,m,v[t-1],r+2*m))
print(f"identity  v_t = r + 2m  checked on {tot} elements, {bad} failures")
for x in ex: print("   ",x)
