"""
Bipartite parity constraint.

The Dynkin diagram is a tree, so 2-colour it c':I->{0,1}.  Every edge changes
the level by 1, so ALL elements of colour i sit at levels of a single parity.
Hence colour i has only ceil/floor of the levels available, and

    v_i  <=  #{ levels of the right parity in [0, h] }

is forced, where h is the height of the poset.  Test whether violating this is
exactly the simple-spectrum failure.
"""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve

def bipart(rs):
    col = {}
    stack = [(rs.I[0], 0)]
    while stack:
        x, c = stack.pop()
        if x in col: continue
        col[x] = c
        for y in rs.adj[x]: stack.append((y, 1 - c))
    return col

for typ, n, ks in [("D",5,[2,3]), ("D",6,[2,3]), ("E",6,[2,3]), ("E",7,[1])]:
    rs = RootSystem(typ, n); col = bipart(rs)
    for k in ks:
        orb = orbit_with_words(rs, rs.fundamental(k))
        if len(orb) > 400: continue
        agree = dis = 0; ex = []
        for mu, word in orb.items():
            if not word: continue
            v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
            r = solve(rs, k, v)
            if not r: continue
            rows = list(r.values())[0]
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            h = max(rows)                      # top level
            # available levels of each parity in [0,h]
            avail = {p: sum(1 for d in range(h+1) if d % 2 == p) for p in (0,1)}
            base = col[k]                      # colour k sits at level 0
            pred = all(v[i-1] <= avail[(col[i]-base) % 2] for i in rs.I)
            if ss == pred: agree += 1
            else:
                dis += 1
                if len(ex) < 3: ex.append((word, v, ss, pred, h))
        print(f"{rs.name()} om{k}: parity bound agrees {agree}, disagrees {dis}")
        for w, v, s, p, h in ex:
            print(f"    w={''.join('s'+str(i) for i in w)} v={v} h={h} simple={s} pred={p}")
