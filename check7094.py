"""Does 'the c-trace attains the value 3' mean the same as 'max c_p = 3'?
If some max-4 element also contains a 3, the theorem's hypothesis is wider
than the 7094 it quotes.
usage: python3 check7094.py
"""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
tot_has3=tot_max3=tot_max4_with3=0
for T_,N in (('E',6),('E',7)):
    rs=RootSystem(T_,N); h3=m3=m4w3=0
    for k in rs.I:
        lam=rs.fundamental(k)
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            cs=[c for _,_,c,_ in trace(rs,w,lam)]
            if 3 in cs: h3+=1
            if max(cs)==3: m3+=1
            if max(cs)>=4 and 3 in cs: m4w3+=1
    print(f"{T_}{N}: contains a 3: {h3};  max c = 3 exactly: {m3};  "
          f"max c >= 4 AND contains a 3: {m4w3}")
    tot_has3+=h3; tot_max3+=m3; tot_max4_with3+=m4w3
print(f"\nE6+E7: contains a 3 = {tot_has3};  max c = 3 = {tot_max3};  "
      f"max>=4 with a 3 = {tot_max4_with3}")
print(f"appendixE quotes 7094.  matches 'max c = 3': {tot_max3==7094}; "
      f"matches 'contains a 3': {tot_has3==7094}")
