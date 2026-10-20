"""Test the two surviving side conditions at EVERY neighbour of t, on every
cached type.  Type D is the control: there, c_t >= 1 must already work.
usage: python3 escand.py d6star.tsv e6star.tsv [e7star.tsv]
"""
import sys, csv
from rootsys import RootSystem
for fn in sys.argv[1:]:
    rows=[]
    for r in csv.DictReader(open(fn),delimiter='\t'):
        rows.append(dict(k=int(r['kappa']),c=[int(x) for x in r['c'].split('|')],
                         sh=r['shape'],ups=set(int(x) for x in r['upsets'].split('|') if x)))
    n=len(rows[0]['c']); rs=RootSystem('E' if fn[0]=='e' else 'D', n)
    t=[i for i in rs.I if len(rs.adj[i])>=3][0]; nbrs=sorted(rs.adj[t])
    leaf=[j for j in nbrs if len(rs.adj[j])==1]
    print(f"\n=== {fn} (node {t}, nbrs {nbrs}, leaf nbrs {leaf}) ===")
    C={'c_t>=1 (type D)': lambda r: r['c'][t-1]>=1,
       'c_t>=2':          lambda r: r['c'][t-1]>=2,
       'c_t>=3':          lambda r: r['c'][t-1]>=3,
       'c_t>=1 & c_leaf>=1': lambda r: r['c'][t-1]>=1 and all(r['c'][j-1]>=1 for j in leaf),
       'c_t>=1 & c_leaf>=0': lambda r: r['c'][t-1]>=1 and all(r['c'][j-1]>=0 for j in leaf)}
    for name,f in C.items():
        good=[r for r in rows if f(r)]
        line=[]
        for j in nbrs:
            bad=sum(1 for r in good if j not in r['ups'])
            line.append(f"j={j}:{bad}")
        ok=all(x.split(':')[1]=='0' for x in line) and good
        print(f"  {name:22} covers {len(good):6}  failures  " + "  ".join(line)
              + ("   *** holds at every neighbour ***" if ok else ""))
