"""E6 omega_3: elements sharing v_4 = 3 that land on OPPOSITE sides.
Look for the invariant that separates them."""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from dfs3 import solve

rs = RootSystem("E", 6); k = 3; lam = rs.fundamental(k); t = 4
S, NS = [], []
for mu, word in sorted(orbit_with_words(rs, lam).items(), key=lambda x: len(x[1])):
    if not word: continue
    tr = trace(rs, word, lam); v = v_from_trace(rs, lam and tr)
    cs = [c for _, _, c, _ in tr]
    r = solve(rs, k, v)
    if not r: continue
    rows = list(r.values())[0]
    ss = all(c <= 1 for x in rows.values() for c in x.values())
    if v[t-1] != 3: continue
    rec = dict(word=word, v=v, cs=cs, cmax=max(cs), ndef=sum(1 for c in cs if c>=2),
               ell=len(word), tot=sum(v), rows=rows)
    (S if ss else NS).append(rec)

print(f"v_4 = 3 elements:  simple {len(S)},  NOT simple {len(NS)}\n")
def show(tag, L):
    print(f"--- {tag}")
    for r in L[:8]:
        print(f"  w={''.join('s'+str(i) for i in r['word']):<30} v={r['v']} "
              f"l={r['ell']} |v|={r['tot']} ndef={r['ndef']} c={r['cs']}")
show("SIMPLE", S); print(); show("NOT SIMPLE", NS)

print("\n--- candidate separators")
for name, f in [
    ("number of defect steps (c>=2)", lambda r: r['ndef']),
    ("length", lambda r: r['ell']),
    ("|v| total", lambda r: r['tot']),
    ("v_2", lambda r: r['v'][1]),
    ("v_3", lambda r: r['v'][2]),
    ("v_5", lambda r: r['v'][4]),
    ("max v_i", lambda r: max(r['v'])),
    ("sum of v over neighbours of t", lambda r: sum(r['v'][j-1] for j in rs.adj[t])),
    ("v_t - max neighbour v", lambda r: r['v'][t-1] - max(r['v'][j-1] for j in rs.adj[t])),
]:
    a = sorted({f(r) for r in S}); b = sorted({f(r) for r in NS})
    sep = "SEPARATES" if not (set(a) & set(b)) else ""
    print(f"  {name:<34} simple={a}  not={b}  {sep}")
