"""
Trace a reduced word acting on lambda, recording how many copies of each
simple root are subtracted at each step.

For w = s_{i_1} ... s_{i_l} acting on lambda (rightmost letter first), set
    mu_0 = lambda,   c_k = <alpha_{i}^vee, mu_{k-1}>,   mu_k = mu_{k-1} - c_k alpha_i.

Then  lambda - w(lambda) = sum_k c_k alpha_{i_k},  so  v_i = sum over steps of color i.

Definition 5.1 (lambda-minuscule) is exactly the condition  c_k = 1 for all k.
When some c_k >= 2 the step contributes c_k elements of the same color to the
poset -- this is where simple spectrum can break.
"""

from rootsys import RootSystem


def trace(rs, word, lam):
    """Returns list of (step, color, c, weight_after), rightmost letter first."""
    mu = tuple(lam)
    out = []
    for step, i in enumerate(reversed(word), start=1):
        c = mu[i - 1]
        mu = tuple(mu[j - 1] - c * rs.A[j - 1][i - 1] for j in rs.I)
        out.append((step, i, c, mu))
    return out


def v_from_trace(rs, tr):
    v = [0] * rs.n
    for _, i, c, _ in tr:
        v[i - 1] += c
    return tuple(v)


def report(label, typ, n, k, word, vexp, poset_desc):
    rs = RootSystem(typ, n)
    lam = rs.fundamental(k)
    tr = trace(rs, word, lam)
    v = v_from_trace(rs, tr)
    print(f"{label}   {rs.name()}  lam=om{k}  w={''.join('s'+str(i) for i in word)}")
    print("    step  letter  c   weight after")
    for step, i, c, mu in tr:
        flag = "   <-- c>1" if c > 1 else ""
        print(f"     {step:2d}     s{i}    {c}   {mu}{flag}")
    print(f"    colour multiset from trace: "
          + " ".join(f"{i}x{v[i-1]}" for i in rs.I if v[i - 1]))
    print(f"    v = {v}   paper v = {vexp}   {'MATCH' if v == vexp else 'DIFFERS'}")
    print(f"    paper poset (bottom-up levels): {poset_desc}")
    print(f"    lambda-minuscule (all c=1): {all(c == 1 for _, _, c, _ in tr)}")
    print()
    return v == vexp


if __name__ == "__main__":
    print("Reduced-word trace: where the extra poset elements come from\n")
    ok = []
    ok.append(report("Ex 5.10", "D", 4, 2, (1, 2, 3, 4, 2), (2, 2, 1, 1),
                     "{2}, {1,3,4}, {2}, {1}"))
    ok.append(report("Ex 5.11", "D", 5, 2, (2, 1, 3, 4, 5, 3, 2), (1, 3, 2, 1, 1),
                     "{2}, {1,3}, {2,4,5}, {3}, {2}"))
    ok.append(report("Ex 5.12", "D", 5, 2, (2, 1, 2, 3, 4, 5, 3, 2), (2, 3, 2, 1, 1),
                     "{2}, {1,3}, {2,4,5}, {3}, {2}, {1}"))
    print("All v reproduced from the word trace." if all(ok) else "MISMATCH")
