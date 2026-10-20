"""Stage 1 only: enumerate the {1,3}-traces with one 3 (no simple-spectrum test).
usage: python3 e3scan.py E7 [cap]"""
import sys
from rootsys import RootSystem
typ = sys.argv[1]; T, N = typ[0], int(typ[1:])
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 40
rs = RootSystem(T, N)
for k in rs.I:
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
    print(f"  lam=om{k}: states {len(seen)}, elements with a 3: {len(hits)}"
          + (f", lengths {min(h[1] for h in hits)}..{max(h[1] for h in hits)}" if hits else ""),
          flush=True)
