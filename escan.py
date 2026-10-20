"""Simple-spectrum verdict for every extremal weight, via soc_sigma (validated in socle.py).
usage: python3 escan.py E6 [kmin] [kmax] [maxlen]"""
import sys, collections, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from socle import soc_sigma

def main():
    typ = sys.argv[1]; T, N = typ[0], int(typ[1:])
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else N
    maxlen = int(sys.argv[4]) if len(sys.argv) > 4 else None
    rs = RootSystem(T, N)
    print(f"{typ}: (k, max_p c_p, simple spectrum) -> count, minimal-length witness", flush=True)
    for k in range(kmin, kmax+1):
        lam = rs.fundamental(k)
        t0 = time.time()
        orb = orbit_with_words(rs, lam, max_len=maxlen)
        tab = collections.Counter(); wit = {}
        for mu, word in sorted(orb.items(), key=lambda x: (len(x[1]), x[1])):
            if not word: continue
            tr = trace(rs, word, lam)
            mc = max(c for _, _, c, _ in tr)
            X, rows, dv = soc_sigma(rs, k, word)
            v = tuple(v_from_trace(rs, tr))
            if dv != v:
                print(f"  !! om{k} w={word} soc dimvec {dv} != v {v}", flush=True); continue
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            key = (mc, ss); tab[key] += 1
            if key not in wit:
                wit[key] = ("".join("s"+str(i) for i in word), v,
                            [c for _,_,c,_ in tr], dict(sorted(rows.items())))
        print(f"  om{k}: |orbit|={len(orb)}  ({time.time()-t0:.1f}s)", flush=True)
        for key in sorted(tab):
            w, v, cs, rows = wit[key]
            print(f"    maxc={key[0]} simple={str(key[1]):5} -> {tab[key]:6}   "
                  f"w={w} v={v} c={cs}", flush=True)
            if key[0] >= 3 and key[1]:
                print(f"        rows={rows}", flush=True)

if __name__ == "__main__":
    main()
