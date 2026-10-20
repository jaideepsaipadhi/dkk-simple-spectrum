"""Does the type-D mechanism at the trivalent node transport to E6?

In type D: the colour-t column of V is the GREEDY TOP TRUNCATION of the
colour-t profile of I_kappa to total v_t (prop:T), and that truncation is
multiplicity-free exactly when v_t <= 2 because the profile is (1,2^{k-1},1).
Test both halves in E6, with type D as the control.
usage: python3 eprof.py [E6|D6 ...]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama, rows_of

def colprof(rs, k, t):
    """multiplicities of colour t in I_kappa, by degree."""
    pr = proj(rs, nakayama(rs)[k])
    out = {}
    for d in range(pr.top + 1):
        out[d] = pr.dim.get((d, t), 0)
    return out

def greedy(prof, m):
    """greedy top truncation to total m: fill from the top degree down."""
    out = {}
    for d in sorted(prof, reverse=True):
        if m <= 0: break
        take = min(prof[d], m); 
        if take: out[d] = take; m -= take
    return out, m

for arg in (sys.argv[1:] or ['D6', 'E6']):
    T, N = arg[0], int(arg[1:])
    rs = RootSystem(T, N); t = [i for i in rs.I if len(rs.adj[i]) >= 3][0]
    KR = [k for k in rs.I if 2 <= k <= N-2] if T == 'D' else list(rs.I)
    agree = tot = 0; mism = None; profs = {}
    for k in KR:
        lam = rs.fundamental(k)
        prof = colprof(rs, k, t); profs[k] = [prof[d] for d in sorted(prof)]
        for mu, w in orbit_with_words(rs, lam).items():
            if not w: continue
            _, rows, dv = soc_sigma(rs, k, w)
            col = {d: row[t] for d, row in rows.items() if t in row}
            g, left = greedy(prof, dv[t-1])
            tot += 1
            if col == g and left == 0: agree += 1
            elif mism is None: mism = (k, w, dv[t-1], col, g)
    print(f"{T}{N} (node {t}): V_t is the greedy top truncation in {agree}/{tot} cases")
    if mism: print(f"    first mismatch: kappa={mism[0]} w={mism[1]} v_t={mism[2]}\n"
                   f"      V_t={mism[3]}\n      greedy={mism[4]}")
    for k in KR: print(f"    colour-{t} profile of I_{k}: {profs[k]}")
