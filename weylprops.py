"""
Classifiers for the three conditions DKK rule out in Section 5.4.2:
fully commutative, minuscule, dominant minuscule.

A Weyl group element is represented by a reduced word (tuple of node labels),
acting on weights via RootSystem.act_word (rightmost letter acts first).
"""

from collections import deque
from rootsys import RootSystem


def length(rs, word):
    """True length of the product, via the sign/inversion-free method:
    reduce the word greedily using the action on the weight lattice."""
    # Use the standard fact: l(w) = number of positive roots sent to negative.
    # Cheaper here: count via the permutation action on all roots.
    roots = all_roots(rs)
    neg = 0
    for r in roots:
        img = rs.act_word(word, r)
        if is_negative(rs, img):
            neg += 1
    return neg


def all_roots(rs):
    """All positive roots, in omega coordinates."""
    simple = [rs.simple_root(i) for i in rs.I]
    seen = set(simple)
    frontier = list(simple)
    while frontier:
        new = []
        for r in frontier:
            for i in rs.I:
                s = rs.act(i, r)
                if s not in seen and not is_negative(rs, s):
                    seen.add(s)
                    new.append(s)
        frontier = new
    return list(seen)


def is_negative(rs, r):
    v = rs.root_coords(r)
    for x in v:
        if x != 0:
            return x < 0
    return False


def is_reduced(rs, word):
    return length(rs, word) == len(word)


def reduced_words(rs, word):
    """All reduced words for the element, via braid + commutation moves."""
    start = tuple(word)
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for k in range(len(u) - 1):
            a, b = u[k], u[k + 1]
            if a == b:
                continue
            if b not in rs.adj[a]:                      # commuting move
                nxt = u[:k] + (b, a) + u[k + 2:]
                if nxt not in seen:
                    seen.add(nxt); q.append(nxt)
            elif k + 2 < len(u) + 0 and k + 2 <= len(u) - 1:
                if u[k + 2] == a:                       # braid move aba -> bab
                    nxt = u[:k] + (b, a, b) + u[k + 3:]
                    if nxt not in seen:
                        seen.add(nxt); q.append(nxt)
    return seen


def is_fully_commutative(rs, word):
    """w is FC iff no reduced word contains s_i s_j s_i with i ~ j."""
    for u in reduced_words(rs, word):
        for k in range(len(u) - 2):
            if u[k] == u[k + 2] and u[k + 1] in rs.adj[u[k]]:
                return False
    return True


def is_lambda_minuscule(rs, word, lam):
    """Definition 5.1: <alpha_{i_k}^vee, s_{i_{k+1}}...s_{i_l}(lam)> = 1 for all k."""
    l = len(word)
    for k in range(l):
        tail = word[k + 1:]
        mu = rs.act_word(tail, lam)
        if mu[word[k] - 1] != 1:
            return False
    return True


def minuscule_witnesses(rs, word, bound=2):
    """Search small weights lam (coeffs in [-bound,bound]) with w lambda-minuscule.

    Returns (any_witness, any_dominant_witness).
    """
    from itertools import product as iproduct
    any_w = None
    any_dom = None
    for coeffs in iproduct(range(-bound, bound + 1), repeat=rs.n):
        lam = tuple(coeffs)
        if is_lambda_minuscule(rs, word, lam):
            if any_w is None:
                any_w = lam
            if all(c >= 0 for c in lam):
                any_dom = lam
                break
    return any_w, any_dom


def report(typ, n, word, k=None):
    rs = RootSystem(typ, n)
    w = tuple(word)
    ell = length(rs, w)
    fc = is_fully_commutative(rs, w) if ell == len(w) else None
    wit, dom = minuscule_witnesses(rs, w)
    s = ''.join('s' + str(i) for i in w)
    print(f"{rs.name():4s} w={s}")
    print(f"     length            : {ell}  (word has {len(w)} letters"
          f"{', REDUCED' if ell == len(w) else ', NOT reduced'})")
    print(f"     fully commutative : {fc}")
    print(f"     minuscule         : {wit is not None}"
          + (f"   witness lam={wit}" if wit else ""))
    print(f"     dominant minusc.  : {dom is not None}"
          + (f"   witness lam={dom}" if dom else ""))
    print()


if __name__ == "__main__":
    print("Checking the classifications DKK assert for Examples 5.10-5.12\n")
    print("Ex 5.10  paper says: minuscule, NOT dominant minuscule")
    report("D", 4, (1, 2, 3, 4, 2))
    print("Ex 5.11  paper says: fully commutative, NOT minuscule")
    report("D", 5, (2, 1, 3, 4, 5, 3, 2))
    print("Ex 5.12  paper says: NOT fully commutative")
    report("D", 5, (2, 1, 2, 3, 4, 5, 3, 2))
