"""Type-E examples with a c-trace entry c_p = 3 (appendix item (3) of Krylov's Sept 14 offer).

c_p = 1 for every p characterises the lambda-minuscule w; c_p = 3 is therefore as far from
the minuscule world as the trace gets in the small ranks. For each such w we report the
c-trace, the dimension vector v, whether the T-action has simple spectrum, and the rows of V
(colours per graded piece, top degree first) in the format of the DKK Hasse diagrams.
usage: python3 etrip.py [E6|E7] [kmin] [kmax] [maxlen]"""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from dfs3 import solve

def rows_str(rows):
    out = []
    for d in sorted(rows, reverse=True):
        cs = []
        for col in sorted(rows[d]):
            cs += [str(col)]*rows[d][col]
        out.append("[" + ",".join(cs) + "]")
    return ",".join(out)

def main():
    typ = sys.argv[1] if len(sys.argv) > 1 else "E6"
    T, N = typ[0], int(typ[1:])
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else N
    maxlen = int(sys.argv[4]) if len(sys.argv) > 4 else None
    rs = RootSystem(T, N)
    print(f"{typ}: elements w with max_p c_p = 3")
    print(f"{'k':>2} {'w':<34} {'v':<22} {'c-trace':<26} {'fc':>3} {'simple':>7}  rows")
    found = 0
    for k in range(kmin, kmax+1):
        lam = rs.fundamental(k)
        for mu, word in sorted(orbit_with_words(rs, lam, max_len=maxlen).items(),
                               key=lambda x: (len(x[1]), x[1])):
            if not word: continue
            tr = trace(rs, word, lam)
            cs = [c for _, _, c, _ in tr]
            if max(cs) != 3: continue
            v = v_from_trace(rs, tr)
            r = solve(rs, k, v)
            if not r: continue
            rows = list(r.values())[0]
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            fc = is_fully_commutative(rs, word)
            w = "".join("s"+str(i) for i in word)
            print(f"{k:>2} {w:<34} {str(v):<22} {str(cs):<26} {str(fc)[0]:>3} "
                  f"{str(ss):>7}  {rows_str(rows)}")
            found += 1
            if found >= 40: print("   ... (truncated at 40)"); return
    if not found: print("   none in range")

if __name__ == "__main__":
    main()
