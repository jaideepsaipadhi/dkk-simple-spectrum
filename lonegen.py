"""Is a lone step at the trivalent node always light, in every type with a
branch node?  Complete check, not a sample.

If t occurs once in a reduced word, the suffix after it is u in W_{I\{t}}, and
    c_p = <alpha_t^vee, u.lam> = <alpha_t^vee,lam> + sum_{j~t} v_j(u),
where lam - u.lam = sum_j v_j(u) alpha_j.  So enumerating the parabolic
W_{I\{t}} settles every lone node step at once -- and that group is small
(72 in E6, 288 in E7, 1440 in E8).
usage: python3 lonegen.py
"""
from rootsys import RootSystem

def parabolic(rs, drop):
    """every element of W_{I\{drop}}, as the set of weights it sends lam to
    is not enough -- we need the group, so BFS on reduced words."""
    gens=[i for i in rs.I if i!=drop]
    seen={tuple(rs.fundamental(k) for k in rs.I):()}   # identity, keyed by action
    # represent u by its action on ALL fundamental weights (faithful)
    frontier=[tuple(rs.fundamental(k) for k in rs.I)]
    while frontier:
        new=[]
        for st in frontier:
            for i in gens:
                im=tuple(rs.act(i,m) for m in st)
                if im not in seen:
                    seen[im]=seen[st]+(i,); new.append(im)
        frontier=new
    return seen

for T,N in (('D',5),('D',6),('D',7),('E',6),('E',7),('E',8)):
    rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    P=parabolic(rs,t)
    worst={}
    for k in rs.I:
        lam=rs.fundamental(k); mx=0; arg=None
        for st,word in P.items():
            ulam=rs.act_word(word,lam)
            c=ulam[t-1]                      # <alpha_t^vee, u.lam>
            if c>mx: mx,arg=c,word
        worst[k]=mx
    ok=all(v<=1 for v in worst.values())
    print(f"{T}{N}: |W_(I-t)| = {len(P):5}, max c over lone node steps by kappa: "
          f"{worst}  -> all <= 1: {ok}")
