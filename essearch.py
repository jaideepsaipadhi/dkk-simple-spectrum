"""Search for the side condition that makes (*o) true at each neighbour of t.
Type D's is c_t >= 1.  Scored on whatever tsv files exist; D6 is the control.
usage: python3 essearch.py d6star.tsv e6star.tsv [e7star.tsv]
"""
import sys, csv, collections
for fn in sys.argv[1:]:
    rows=[]
    for r in csv.DictReader(open(fn),delimiter='\t'):
        rows.append(dict(k=int(r['kappa']),
                         c=[int(x) for x in r['c'].split('|')],
                         v=[int(x) for x in r['v'].split('|')],
                         sh=r['shape'],
                         ups=set(int(x) for x in r['upsets'].split('|') if x)))
    n=len(rows[0]['c']); t=4 if n in (6,7,8) else n-2
    # recover t as the branch node by degree, from the diagram
    from rootsys import RootSystem
    T_='E' if n in (6,7,8) and fn[1]=='6' or fn[0]=='e' else 'D'
    rs=RootSystem('E' if fn[0]=='e' else 'D', n)
    t=[i for i in rs.I if len(rs.adj[i])>=3][0]; nbrs=sorted(rs.adj[t])
    print(f"\n=== {fn}  (node {t}, neighbours {nbrs}, {len(rows)} weights) ===")
    for j in nbrs:
        fails=[r for r in rows if j not in r['ups']]
        dist=collections.Counter(r['c'][t-1] for r in fails)
        print(f"  j={j}: {len(fails)} non-up-sets; c_t distribution among them: "
              f"{dict(sorted(dist.items()))}")
        # minimal threshold on c_t that excludes every failure
        if fails:
            mx=max(dist)
            print(f"         -> 'c_t >= {mx+1}' is the weakest c_t-threshold with no failure"
                  f"   (type D's is c_t >= 1)")
        # does a condition at j itself do better?
        for name,cond in [('c_j >= 1', lambda r: r['c'][j-1]>=1),
                          ('c_j >= 0', lambda r: r['c'][j-1]>=0),
                          ('c_t >= 1 and c_j >= 0', lambda r: r['c'][t-1]>=1 and r['c'][j-1]>=0),
                          ('c_t >= 1 and c_j >= 1', lambda r: r['c'][t-1]>=1 and r['c'][j-1]>=1),
                          ('c_t + c_j >= 1', lambda r: r['c'][t-1]+r['c'][j-1]>=1),
                          ('c_t >= 1 and shape full', lambda r: r['c'][t-1]>=1 and r['sh']=='full')]:
            bad=[r for r in rows if cond(r) and j not in r['ups']]
            good=[r for r in rows if cond(r)]
            if not bad and good:
                print(f"         *** '{name}' => up-set at {j}:  {len(good)} cases, 0 failures")
