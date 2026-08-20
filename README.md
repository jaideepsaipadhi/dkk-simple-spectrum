# Extremal torus fixed points and simple spectrum

Code and manuscript answering the open problem of §5.4.2 of Dinkins–Karpov–Krylov,
*The quantum Hikita conjecture via quasimaps*, [arXiv:2608.16746](https://arxiv.org/abs/2608.16746):

> "We consider it an interesting open problem to find a combinatorial
> characterization of such *w*."

## The result

Let Γ be a simply-laced Dynkin diagram, λ = ω_κ a fundamental weight, w ∈ W and
μ = wλ. Then M(v, w) is a point with a unique T-fixed representation V. Let Π be
the preprojective algebra of Γ and I_κ the injective envelope of the simple S_κ
(so I_κ = Πe_{ν⁻¹(κ)} for the Nakayama permutation ν).

**Main theorem.** Read a reduced word for w from the right, recording at each
step c = ⟨α_i^∨, μ_current⟩; let σ be the resulting sequence of vertices with
multiplicities. Then

    V  ≅  soc_σ(I_κ)

as graded modules, where soc_σ is the Geiss–Leclerc–Schröer iterated socle.
Simple spectrum ⟺ no two elements of the resulting coloured graded set share a
colour and a degree — decidable in polynomial time from any reduced word.

The proof is uniform in type. It rests on an identity obtained from
Crawley-Boevey's Ext¹ formula,

    dim soc_i(I_κ/X) = ⟨α_i^∨, μ'⟩ + dim top_i(X),

together with a squeeze against dim I_κ = λ − w₀λ.

**Type D closed form.** For λ = ω_k with k ≤ n−2, write μ = Σ_{i∈S} η_i e_i in
Bourbaki coordinates, m = #{minus signs among the first n−2 coordinates},
r = |S ∩ {n−1, n}|. Then simple spectrum ⟺ **r + 2m ≤ 2**. Unconditional for
the spin weights; machine-checked ingredients for k ≥ 2 — see the table below.

**Reduction to one vertex.** Only the trivalent node t needs checking: if V_t is
multiplicity-free then so is V. Proved in both types, via inward injectivity on
I_κ — exhaustively in type E, and by a rank-uniform argument in type D.

## Quick start

    python3 proof.py        # every step of the proof of the main theorem
    python3 verify.py       # DKK's printed posets, criterion vs solver, type D
    python3 sweepall.py     # aggregate criterion-vs-solver sweep (2179/2179)
    python3 inwardlemma.py  # the inward-injectivity lemma
    python3 master.py       # THE open conjecture, over all graded submodules
    python3 filltop.py      # the unconditional reductions (a) and (b)
    python3 claimG.py       # reduction (c), link by link
    python3 inwardproof.py  # inputs to the type-D proof of Lemma I
    python3 bigrade.py      # the bidegree is a function of (vertex, degree)

Each script prints a pass/fail line and exits 0/1. No dependencies beyond the
standard library. Linear algebra is over F_p; set `DKKP` to change the
characteristic (5, 7, 11, 13 all tested) and `DKKCAP` to raise the orbit cutoff
in `sweepall.py`.

## The manuscript

`paper.tex` / `paper.pdf` (11 pages).

## What is proved and what is not

| statement | status |
|---|---|
| Main theorem, V ≅ soc_σ(I_κ), all ADE | **proved** |
| Crawley-Boevey socle-step identity | **proved** |
| T-weight grading = path-length grading | **proved** |
| dim I_κ = λ − w₀λ | classical; verified 47/47 |
| minuscule λ ⟹ always simple | **proved** |
| checking the trivalent node suffices, type E | **proved** (exhaustive finite check) |
| checking the trivalent node suffices, type D | **proved** (uniformly in the rank) |
| type D closed form r + 2m ≤ 2, spin weights | **proved** |
| type D closed form, k ≥ 2, leaf-symmetric v | **proved** (via the τ-twist) |
| type D closed form, k ≥ 2, otherwise | rests on the staircase conjecture |

**The one open point.** Call a graded submodule M ⊆ I_κ a **staircase** if its
colour-t multiplicities are the greedy top truncation of those of I_κ (zero
below one degree, proper at it, full above) and at each leaf M occupies a final
segment of degrees. Everything left open is:

> **Conjecture.** If dim M is w-extremal then M is a staircase. Equivalently, by
> Baumann–Kamnitzer–Tingley rigidity, a non-staircase submodule has
> δ(dim M) = 2v_κ − vᵀCv > 0 and so is never the unique submodule of its
> dimension vector.

It gives Proposition T at once: an extremal v carries a unique submodule, namely
V = soc_σ(I_κ), and the conjecture makes it a staircase. Verified on all 1124
graded submodules of I_κ over D₄–D₇ — exactly 78 non-staircases, every one with
δ = 2 (`master.py`).

**Proved: the conjecture holds whenever v_{n−1} = v_n.** If M fails at t, every
nonzero x ∈ M_t[d] must be special (annihilated by an arrow to a leaf), else its
double paths would span (I_κ)_t[d+2]. The special vectors are two *lines*, not a
subspace, so M_t[d] is one of them. The leaf-exchanging automorphism τ fixes I_κ
and exchanges the two special lines, so τ(M) ≠ M; when v_{n−1} = v_n it
preserves dim M, and two submodules of one dimension vector contradict rigidity.

**Three unconditional reductions** (running the induction along socle steps,
D = colour-t threshold of the partial module X):

- (a) a colour-t step contributes nothing below D−2;
- (b) if it contributes at D−2 then dim X_t[D] = 1 and that contribution is
  exactly one special line, so one further dimension at D suffices;
- (c) assuming the leaf clause, both leaves are automatically full at D+1:
  c = μ'_t − μ'_{n−1} ≥ 1 with μ'_t ≤ 1 forces μ'_{n−1} ≤ 0, and then
  v'_{n−1} + v'_n ∈ {2a, 2a+1} with |v'_{n−1} − v'_n| ≤ 1 gives
  min(v'_{n−1}, v'_n) ≥ a.

(a) and (b) follow from the simple socle plus the Lemma I propagation.

**Two facts that delimit the search.** The conjecture genuinely needs
extremality — the staircase property fails for 78 of those 1124 submodules. And
the obvious mechanism is unavailable: when X_t[D] is a line it is *always* one
of the two special lines (527/527), so one leaf arrow kills it and that leaf's
fullness at D+1 comes from elsewhere in X.

Two things that were open in earlier drafts and are now closed:

- *Checking only the trivalent node* is proved in both types via inward
  injectivity on I_κ. In type D there is a rank-uniform proof (propagate the
  kernel element outward along the chain using the relation and the simple
  socle; at vertex 1 it lands in the socle, forcing κ = 1, which Lemma P
  excludes). The obvious argument — containment of multiplicity-2 windows —
  genuinely fails, and not only in type E: in D₆ with κ=3 vertex 3 carries
  multiplicity 2 at degree 4 while the trivalent node carries it at 3 and 5.
- *The grading convention* is settled: Γ is a tree, so on Πe_j the bidegree
  (#a, #a*) is determined by the vertex and the path-length degree, and any
  torus acting through the bigrading sees a weight that is an injective affine
  function of degree. So the criterion is convention-independent.

## Problem 2 (separate)

`krylov-wang-conj16/` holds an independent computational attack on Conjecture
1.6 of Krylov–Wang, [arXiv:2608.03314](https://arxiv.org/abs/2608.03314) — a
different problem, sharing no code with the above. See its own README.

## Core files

| file | contents |
|---|---|
| `rootsys.py` | ADE root systems, Weyl action, orbits with reduced words |
| `wordtrace.py` | the c-trace |
| `proj2.py` | Πe_k built without walk enumeration |
| `gls.py` | the GLS iterated socle |
| `dfs3.py` | independent graded-submodule solver (the ground truth) |
| `homological.py` | the socle-step identity |
| `squeeze.py` | dim I_κ = λ − w₀λ |
| `epsilon.py` | ε-coordinates, v_t = r + 2m |
| `inwardlemma.py` | inward injectivity (Lemma I) |
| `inwardproof.py` | the inputs to the rank-uniform type-D proof of Lemma I |
| `bigrade.py` | Γ is a tree, so the bidegree is determined by (vertex, degree) |
| `violators.py` | δ, |supp μ| and leaf symmetry of the upward-closure violators |
| `filltop.py` | the socle-step induction reducing Proposition T to statement (b) |
| `epscrit.py` | the ε-coordinate signature of the critical steps |
| `claimG.py` | Claim G and the proved leaf case, link by link |
| `upclosure.py` | Lemma U, Lemma S′, the defect δ |

`RESULT.md`, `THEOREM.md`, `FINDINGS.md` are the working research logs, kept for
the record; they contain superseded formulations and are not the reference.
