"""A PROVED bound on the c-trace under simple spectrum, valid in every type
including E8, where the Weyl orbits cannot be enumerated.

  v_i = sum of c_p over the occurrences of colour i        (eq:vsumc, proved)
  simple spectrum => the colour-i column of V is multiplicity-free
                  => v_i <= L_i, the number of nonzero degrees of the
                     colour-i column of I_kappa
  hence            c_p <= v_{i_p} <= L_{i_p}  for every p.

So max_i L_i bounds every c-trace value of every simple-spectrum element.
usage: python3 cbound.py
"""
from rootsys import RootSystem
from socle import proj, nakayama
for T,N in (('D',5),('D',6),('D',7),('E',6),('E',7),('E',8)):
    rs=RootSystem(T,N)
    print(f"\n{T}{N}:")
    worst=0
    for k in rs.I:
        pr=proj(rs,nakayama(rs)[k])
        L={i:sum(1 for (d,j) in pr.dim if j==i and pr.dim[(d,i)]>0) for i in rs.I}
        m=max(L.values()); worst=max(worst,m)
        print(f"   kappa={k}: column lengths L_i = {[L[i] for i in rs.I]}   bound max_i L_i = {m}")
    print(f"   => c_p <= {worst} for every simple-spectrum element of {T}{N}")
