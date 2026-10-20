"""Stage 2: test the E7 {1,3}-trace candidates for simple spectrum, one at a time."""
import sys, time
from rootsys import RootSystem
from dfs3 import solve
typ = sys.argv[1]; T, N = typ[0], int(typ[1:])
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 40
rs = RootSystem(T, N)
KS = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv)>3 else [])] or rs.I
for k in KS:
    lam = rs.fundamental(k); hits = []; seen = set()
    st = [(tuple(lam), 0, 0, ())]
    while st:
        mu, n3, L, word = st.pop()
        key = (mu, n3)
        if key in seen: continue
        seen.add(key)
        if n3: hits.append((mu, L, word))
        if L >= cap: continue
        for i in rs.I:
            c = mu[i-1]
            if c == 1: st.append((rs.act(i, mu), n3, L+1, word+(i,)))
            elif c == 3 and n3 == 0: st.append((rs.act(i, mu), 1, L+1, word+(i,)))
    for mu, L, word in sorted(hits, key=lambda h: h[1]):
        v = tuple(int(x) for x in rs.root_coords(tuple(lam[j]-mu[j] for j in range(N))))
        t0 = time.time()
        r = solve(rs, k, v)
        dt = time.time() - t0
        w = "".join("s"+str(i) for i in reversed(word))
        if not r:
            print(f"  om{k} L={L} v={v}  NO SOLUTION  ({dt:.1f}s)  w={w}", flush=True); continue
        rows = list(r.values())[0]
        ss = all(c <= 1 for x in rows.values() for c in x.values())
        rs_ = ",".join("[" + ",".join(str(c) for c in sorted(rows[d]) for _ in range(rows[d][c]))
                       + "]" for d in sorted(rows, reverse=True))
        print(f"  om{k} L={L} mu={mu}\n      v={v} simple={ss} ({dt:.1f}s)\n"
              f"      w={w}\n      rows={rs_}", flush=True)
