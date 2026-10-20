"""Three structural checks behind a tighter statement of thm:typeD.

(1) v_i = sum of c_p over occurrences of colour i  (sigma is the c-weighted
    sequence, so this should be an identity, not a conjecture).
(2) simple spectrum <=> v_t <= 2, t the trivalent node -- the basis-free form.
(3) c_p >= 2 at some occurrence of t  ==>  t occurs at least twice,
    hence v_t >= 3.  This is the bridge: it makes the trivalent-node
    necessary condition a COROLLARY of (2) rather than a separate fact.
"""
import collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma

for T,N in (('D',4),('D',5),('D',6),('D',7),('E',6)):
    rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    bad1=bad2=bad3=tot=0; ex={}
    for k in rs.I:
        lam=rs.fundamental(k)
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            tr=trace(rs,w,lam); tot+=1
            _,rows,dv=soc_sigma(rs,k,w)
            simple=all(c<=1 for x in rows.values() for c in x.values())
            # (1)
            s=collections.Counter()
            for _,i,c,_ in tr: s[i]+=c
            if any(s[i]!=dv[i-1] for i in rs.I):
                bad1+=1; ex.setdefault(1,(k,w,dict(s),dv))
            # (2) only claimed in type D
            if T=='D' and simple != (dv[t-1]<=2):
                bad2+=1; ex.setdefault(2,(k,w,dv,simple))
            # (3)
            occ=[c for _,i,c,_ in tr if i==t]
            if any(c>=2 for c in occ):
                if len(occ)<2 or dv[t-1]<3:
                    bad3+=1; ex.setdefault(3,(k,w,occ,dv))
    print(f"{T}{N} (node {t}), {tot} extremal weights: "
          f"(1) v_i=sum c_p fails {bad1}; "
          f"(2) simple<=>v_t<=2 fails {bad2}{' [n/a]' if T!='D' else ''}; "
          f"(3) bridge fails {bad3}")
    for j in sorted(ex): print(f"     counterexample ({j}): {ex[j]}")
