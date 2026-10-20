"""Type E, elements w whose c-trace attains 3: does simple spectrum ever survive?
Reports, per (k, max c), counts of simple / not simple, with a minimal-length witness of each.
usage: python3 etab.py [E6|E7] [kmin] [kmax] [maxlen]"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from dfs3 import solve

def main():
    typ = sys.argv[1] if len(sys.argv) > 1 else "E6"
    T, N = typ[0], int(typ[1:])
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else N
    maxlen = int(sys.argv[4]) if len(sys.argv) > 4 else None
    rs = RootSystem(T, N)
    tab = collections.Counter(); wit = {}
    for k in range(kmin, kmax+1):
        lam = rs.fundamental(k)
        orb = orbit_with_words(rs, lam, max_len=maxlen)
        for mu, word in sorted(orb.items(), key=lambda x: (len(x[1]), x[1])):
            if not word: continue
            tr = trace(rs, word, lam)
            mc = max(c for _, _, c, _ in tr)
            v = v_from_trace(rs, tr)
            r = solve(rs, k, v)
            if not r: continue
            rows = list(r.values())[0]
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            key = (k, mc, ss)
            tab[key] += 1
            if key not in wit:
                wit[key] = ("".join("s"+str(i) for i in word), tuple(v),
                            [c for _,_,c,_ in tr], is_fully_commutative(rs, word))
    print(f"{typ}: (k, max_p c_p, simple spectrum) -> count, minimal-length witness")
    for key in sorted(tab):
        w, v, cs, fc = wit[key]
        print(f"  k={key[0]} maxc={key[1]} simple={str(key[2]):5} -> {tab[key]:6}   "
              f"w={w} v={v} c={cs} fc={fc}")

if __name__ == "__main__":
    main()
