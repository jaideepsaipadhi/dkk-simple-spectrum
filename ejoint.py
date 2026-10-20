"""Predict 'strand mode' from (kappa, v), scored on E6 AND E7 JOINTLY.
Features are defined by diagram position (arm neighbour / arm far end / short
arm / node), not by vertex number, so they mean the same thing in both.
usage: python3 ejoint.py
"""
import csv, collections, itertools
D={}
for tag,fn,arms in (('E6','e6shape.tsv',{'L1':[3,1],'S':[2],'L2':[5,6]}),
                    ('E7','e7shape.tsv',{'L1':[3,1],'S':[2],'L2':[5,6,7]})):
    rows=[]
    for r in csv.DictReader(open(fn),delimiter='\t'):
        v=[int(x) for x in r['v'].split('|')]
        rows.append(dict(ka=int(r['kappa']),v=v,shape=r['shape'],
                         simple=int(r['simple']),vt=int(r['vt'])))
    D[tag]=(rows,arms)
# only the band matters: shape 'strand' never occurs with v_t small
for tag,(rows,arms) in D.items():
    sh=collections.Counter((r['shape'],r['vt']-0) for r in rows if r['shape']=='strand')
    print(f"{tag}: strand cases by v_t -> {dict(collections.Counter(r['vt'] for r in rows if r['shape']=='strand'))}")

def feats(tag,r):
    rows,arms=D[tag]
    v=r['v']; k=r['ka']
    n1,f1=arms['L1'][0],arms['L1'][-1]
    n2,f2=arms['L2'][0],arms['L2'][-1]
    s=arms['S'][0]
    armof={}
    for nm,a in arms.items():
        for i in a: armof[i]=nm
    ak=armof.get(k,'NODE')
    return dict(nb1=v[n1-1],nb2=v[n2-1],f1=v[f1-1],f2=v[f2-1],sh=v[s-1],
                arm=ak,ka=k,vt=r['vt'])

CAND={
 'nb parities both even':      lambda F: F['nb1']%2==0 and F['nb2']%2==0,
 'far-end parity = not-own-arm': lambda F: (F['f1']%2==(0 if F['arm']=='L1' else 1)
                                        and F['f2']%2==(0 if F['arm']=='L2' else 1)),
 'both of the above':          lambda F: (F['nb1']%2==0 and F['nb2']%2==0
                                        and F['f1']%2==(0 if F['arm']=='L1' else 1)
                                        and F['f2']%2==(0 if F['arm']=='L2' else 1)),
 'short-arm even':             lambda F: F['sh']%2==0,
 'nb1+nb2 even':               lambda F: (F['nb1']+F['nb2'])%2==0,
 'f1+f2 parity vs arm':        lambda F: (F['f1']+F['f2'])%2==(1 if F['arm']=='NODE' else 0),
 'vt even':                    lambda F: F['vt']%2==0,
 'vt odd':                     lambda F: F['vt']%2==1,
}
print()
for name,f in CAND.items():
    line=[]
    for tag in ('E6','E7'):
        rows,_=D[tag]
        tp=fp=fn=0
        for r in rows:
            F=feats(tag,r); p=f(F); s=(r['shape']=='strand')
            if p and s: tp+=1
            elif p and not s: fp+=1
            elif s: fn+=1
        line.append(f"{tag}: tp={tp} fp={fp} fn={fn}")
    print(f"{name:34} " + "   ".join(line))
