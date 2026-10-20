"""Is the column shape a function of LOCAL data at the node?

A submodule can only take a line at a multiplicity->=2 node slot if the arms
deliver one dimension there, so the shape should be forced by v at t and its
neighbours.  Test whether shape is constant on each group.
usage: python3 elocal.py
"""
import csv, collections
NB={'E6':[2,3,5],'E7':[2,3,5]}
for tag,fn in (('E6','e6shape.tsv'),('E7','e7shape.tsv')):
    rows=[]
    for r in csv.DictReader(open(fn),delimiter='\t'):
        v=[int(x) for x in r['v'].split('|')]
        rows.append((int(r['kappa']),v,r['shape'],int(r['simple'])))
    print(f"\n{tag}: {len(rows)} weights")
    for name,key in [
      ('(kappa, v_t, v at nbrs 2,3,5)', lambda k,v:(k,v[3],v[1],v[2],v[4])),
      ('(v_t, v at nbrs 2,3,5)',        lambda k,v:(v[3],v[1],v[2],v[4])),
      ('(kappa, v_t, v_3, v_5)',        lambda k,v:(k,v[3],v[2],v[4])),
      ('(kappa, v at nbrs only)',       lambda k,v:(k,v[1],v[2],v[4])),
      ('full v (control: must be 1-1)', lambda k,v:(k,tuple(v))),
    ]:
        g=collections.defaultdict(set)
        for k,v,sh,s in rows: g[key(k,v)].add(sh)
        mixed=sum(1 for s in g.values() if len(s)>1)
        print(f"   {name:34} {len(g):6} groups, {mixed:5} with >1 shape"
              + ("   *** SHAPE IS DETERMINED ***" if mixed==0 else ""))
