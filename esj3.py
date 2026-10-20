"""What side condition makes (*o) true at j=3 in E6?  The 16 failures with
c_t >= 1 all have c_t = 1; 12 of them have FULL column shape, so
'c_t>=1 and full' does not do it.  Wider search.
usage: python3 esj3.py e6star.tsv 3
"""
import sys, csv, collections, itertools
fn=sys.argv[1]; J=int(sys.argv[2])
rows=[]
for r in csv.DictReader(open(fn),delimiter='\t'):
    rows.append(dict(k=int(r['kappa']),c=[int(x) for x in r['c'].split('|')],
                     v=[int(x) for x in r['v'].split('|')],sh=r['shape'],
                     ups=set(int(x) for x in r['upsets'].split('|') if x)))
t=4
F=[r for r in rows if J not in r['ups']]
P=[r for r in F if r['c'][t-1]>=1]
print(f"{fn}: j={J}, {len(F)} non-up-sets, {len(P)} with c_t>=1")
print("the c_t>=1 failures:")
for r in P: print(f"   kappa={r['k']} c={r['c']} v={r['v']} shape={r['sh']}")
# candidate side conditions, each must give ZERO failures at j=J
C={
 'c_t>=1': lambda r: r['c'][t-1]>=1,
 'c_t>=1 & c_J>=0': lambda r: r['c'][t-1]>=1 and r['c'][J-1]>=0,
 'c_t>=1 & c_J>=1': lambda r: r['c'][t-1]>=1 and r['c'][J-1]>=1,
 'c_t>=1 & c_2>=0': lambda r: r['c'][t-1]>=1 and r['c'][1]>=0,
 'c_t>=1 & c_2>=1': lambda r: r['c'][t-1]>=1 and r['c'][1]>=1,
 'c_t>=1 & c_other>=0': lambda r: r['c'][t-1]>=1 and r['c'][(8-J)-1]>=0,
 'c_t>=1 & c_J+c_other>=0': lambda r: r['c'][t-1]>=1 and r['c'][J-1]+r['c'][(8-J)-1]>=0,
 'c_t>=2': lambda r: r['c'][t-1]>=2,
 'c_t>=1 & all c_nbr>=0': lambda r: r['c'][t-1]>=1 and all(r['c'][i-1]>=0 for i in (2,3,5)),
 'c_t>=1 & v_J even': lambda r: r['c'][t-1]>=1 and r['v'][J-1]%2==0,
 'c_t>=1 & kappa!=4': lambda r: r['c'][t-1]>=1 and r['k']!=4,
 'c_t>=1 & c_J>=0 & c_2>=0': lambda r: r['c'][t-1]>=1 and r['c'][J-1]>=0 and r['c'][1]>=0,
}
print()
for name,f in C.items():
    good=[r for r in rows if f(r)]
    bad=[r for r in good if J not in r['ups']]
    flag='*** 0 failures ***' if not bad and good else f'{len(bad)} failures'
    print(f"  {name:28} covers {len(good):5}  {flag}")
