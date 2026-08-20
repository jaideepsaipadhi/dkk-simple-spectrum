"""
Root systems / Weyl groups for simply-laced (ADE) Dynkin diagrams.

Conventions follow Dinkins-Karpov-Krylov, arXiv:2608.16746:
  D_n : chain 1 - 2 - ... - (n-2), with (n-2) also joined to (n-1) and n.
        So vertex n-2 is the trivalent node.
        (D4: vertex 2 trivalent, edges 1-2,2-3,2-4  -- matches Sec 5.4.3)
        (D5: edges 1-2-3, 3-4, 3-5                  -- matches Sec 5.4.4)

Weights are stored in the FUNDAMENTAL WEIGHT basis:
    a weight mu is the tuple ( <alpha_i^vee, mu> )_{i in I}.
In this basis the simple root alpha_i has coordinates = column i of the
Cartan matrix A, and the simple reflection acts by

    s_i(mu)_j = mu_j - mu_i * A[j][i].
"""

from itertools import product
from sympy import Matrix, Rational


# ---------------------------------------------------------------- diagrams

def edges_A(n):
    return [(i, i + 1) for i in range(1, n)]


def edges_D(n):
    assert n >= 4
    e = [(i, i + 1) for i in range(1, n - 2)]      # 1-2-...-(n-2)
    e += [(n - 2, n - 1), (n - 2, n)]              # trivalent at n-2
    return e


def edges_E(n):
    # Bourbaki E_n: chain 1-3-4-5-...-n with 2 attached to 4
    assert n in (6, 7, 8)
    e = [(1, 3), (3, 4), (2, 4)]
    e += [(i, i + 1) for i in range(4, n)]
    return e


def diagram(typ, n):
    return {"A": edges_A, "D": edges_D, "E": edges_E}[typ](n)


class RootSystem:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        self.I = list(range(1, n + 1))
        self.edges = diagram(typ, n)
        self.adj = {i: set() for i in self.I}
        for a, b in self.edges:
            self.adj[a].add(b)
            self.adj[b].add(a)
        # Cartan matrix, 0-indexed internally
        A = [[0] * n for _ in range(n)]
        for i in self.I:
            A[i - 1][i - 1] = 2
            for j in self.adj[i]:
                A[i - 1][j - 1] = -1
        self.A = A
        self.Amat = Matrix(A)
        self.Ainv = self.Amat.inv()

    def name(self):
        return f"{self.typ}{self.n}"

    # -- weights in fundamental-weight basis -------------------------------

    def fundamental(self, k):
        return tuple(1 if i == k else 0 for i in self.I)

    def simple_root(self, k):
        """alpha_k in the fundamental weight basis = column k of A."""
        return tuple(self.A[i - 1][k - 1] for i in self.I)

    def act(self, i, mu):
        """s_i applied to a weight given in fundamental weight coordinates."""
        c = mu[i - 1]
        if c == 0:
            return mu
        return tuple(mu[j - 1] - c * self.A[j - 1][i - 1] for j in self.I)

    def act_word(self, word, mu):
        """Apply s_{w1} s_{w2} ... s_{wl} to mu (rightmost acts first)."""
        for i in reversed(word):
            mu = self.act(i, mu)
        return mu

    def root_coords(self, mu):
        """Write mu (given in omega-basis) as sum v_i alpha_i; returns tuple of v_i.

        Since alpha has omega-coords A*v, we solve A v = mu.
        """
        v = self.Ainv * Matrix(list(mu))
        return tuple(v[i, 0] for i in range(self.n))

    def pairing(self, mu, nu):
        """Symmetric form (mu,nu), normalised so (alpha_i,alpha_i)=2.

        In omega coordinates (mu,nu) = mu^T A^{-1} nu.
        """
        return (Matrix([list(mu)]) * self.Ainv * Matrix(list(nu)))[0, 0]


# ---------------------------------------------------------------- Weyl orbit

def orbit_with_words(rs, lam, max_len=None):
    """BFS over the Weyl orbit of lam.

    Returns dict: weight -> a shortest word w (as tuple of node labels)
    with w(lam) = weight.
    """
    seen = {tuple(lam): ()}
    frontier = [tuple(lam)]
    while frontier:
        new = []
        for mu in frontier:
            if max_len is not None and len(seen[mu]) >= max_len:
                continue
            for i in rs.I:
                nu = rs.act(i, mu)
                if nu not in seen:
                    seen[nu] = (i,) + seen[mu]
                    new.append(nu)
        frontier = new
    return seen


# ---------------------------------------------------------------- self-test

def _check(typ, n, k, word, exp_mu_desc, exp_v):
    """Reproduce a worked example from the paper."""
    rs = RootSystem(typ, n)
    lam = rs.fundamental(k)
    mu = rs.act_word(word, lam)
    v = rs.root_coords(tuple(lam[i] - mu[i] for i in range(n)))
    dimX = rs.pairing(lam, lam) - rs.pairing(mu, mu)
    ok = tuple(v) == tuple(Rational(x) for x in exp_v)
    print(f"  {rs.name():4s} lam=om{k}  w={''.join('s'+str(i) for i in word)}")
    print(f"       mu (omega coords) = {mu}   [paper: {exp_mu_desc}]")
    print(f"       v = {tuple(int(x) for x in v)}   expected {exp_v}   "
          f"{'OK' if ok else 'MISMATCH'}")
    print(f"       dim X = {dimX}  (expect 0 since mu in W.lam)")
    return ok


if __name__ == "__main__":
    print("Reproducing Examples 5.10-5.12 of arXiv:2608.16746\n")
    r = []
    # Ex 5.10: D4, v=(2,2,1,1), w=(0,1,0,0), lam=om2, mu=-alpha1, w=s1s2s3s4s2
    r.append(_check("D", 4, 2, (1, 2, 3, 4, 2), "-alpha_1", (2, 2, 1, 1)))
    # Ex 5.11: D5, v=(1,3,2,1,1), lam=om2, mu=-alpha2, w=s2s1s3s4s5s3s2
    r.append(_check("D", 5, 2, (2, 1, 3, 4, 5, 3, 2), "-alpha_2", (1, 3, 2, 1, 1)))
    # Ex 5.12: D5, v=(2,3,2,1,1), lam=om2, mu=-alpha1-alpha2, w=s2s1s2s3s4s5s3s2
    r.append(_check("D", 5, 2, (2, 1, 2, 3, 4, 5, 3, 2), "-alpha_1-alpha_2",
                    (2, 3, 2, 1, 1)))
    print("\nAll three reproduced." if all(r) else "\nSOME MISMATCHES")
