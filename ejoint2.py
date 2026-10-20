"""Same joint search, scored only inside the band v_t = T(kappa)+1, which is
where strand mode actually lives in both E6 and E7.
usage: python3 ejoint2.py
"""
import csv, collections
TH={'E6':{1:4,2:3,3:2,4:2,5:2,6:4},'E7':{1:4,2:3,3:2,4:2,5:2,6:2,7:6}}
ARMS={'E6':{'L1':[3,1],'S':[2],'L2':[5,6]},'E7':{'L1':[3,1],'S':[2],'L2':[5,6,7]}}
D={}
for tag,fn in (('E6','e6shape.tsv'),('E7','e7shape.tsv')):
    rows=[]
    for r in csv.DictReader(open(fn),delimiter='\t'):
        v=[int(x) for x in r['v'].split('|')]; k=int(r['kappa'])
        if int(r['vt'])!=TH[tag][k]+1: continue
        rows.append(dict(ka=k,v=v,strand=(r['shape']=='strand'),vt=int(r['vt'])))
    D[tag]=rows
    print(f"{tag} band: {len(rows)} weights, {sum(r['strand'] for r in rows)} strand")
def F(tag,r):
    A=ARMS[tag]; v=r['v']; k=r['ka']
    armof={i:nm for nm,a in A.items() for i in a}
    return dict(nb1=v[A['L1'][0]-1],nb2=v[A['L2'][0]-1],
                f1=v[A['L1'][-1]-1],f2=v[A['L2'][-1]-1],
                sh=v[A['S'][0]-1],arm=armof.get(k,'NODE'),vt=r['vt'],v=v,ka=k)
C={
 'nb both even':  lambda f: f['nb1']%2==0 and f['nb2']%2==0,
 'far-end parity': lambda f:(f['f1']%2==(0 if f['arm']=='L1' else 1)
                         and f['f2']%2==(0 if f['arm']=='L2' else 1)),
 'nb even AND far-end parity': lambda f:(f['nb1']%2==0 and f['nb2']%2==0
                         and f['f1']%2==(0 if f['arm']=='L1' else 1)
                         and f['f2']%2==(0 if f['arm']=='L2' else 1)),
 'nb even AND short-arm par':  lambda f:(f['nb1']%2==0 and f['nb2']%2==0
                         and f['sh']%2==(1 if f['arm']=='S' else 0)),
 'all three':     lambda f:(f['nb1']%2==0 and f['nb2']%2==0
                         and f['f1']%2==(0 if f['arm']=='L1' else 1)
                         and f['f2']%2==(0 if f['arm']=='L2' else 1)
                         and f['sh']%2==(1 if f['arm']=='S' else 0)),
}
print()
for name,fn in C.items():
    out=[]
    for tag in ('E6','E7'):
        tp=fp=fn_=0
        for r in D[tag]:
            p=fn(F(tag,r))
            if p and r['strand']: tp+=1
            elif p: fp+=1
            elif r['strand']: fn_+=1
        out.append(f"{tag}: tp={tp:3} fp={fp:4} fn={fn_:3}")
    print(f"{name:30} " + "   ".join(out))
