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
k ≤ 2 and for the spin weights; see the caveat below for k ≥ 3.

**Reduction to one vertex.** Only the trivalent node t needs checking: if V_t is
multiplicity-free then so is V. Proved in type D by window containment and in
type E by an inward-injectivity lemma whose verification there is exhaustive.

## Quick start

    python3 proof.py        # every step of the proof of the main theorem
    python3 verify.py       # DKK's printed posets, criterion vs solver, type D
    python3 sweepall.py     # aggregate criterion-vs-solver sweep (2179/2179)
    python3 inwardlemma.py  # the inward-injectivity lemma
    python3 upclosure.py    # state of the k >= 3 truncation argument

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
| dim I_κ = λ − w₀λ | classical; verified 47/47 |
| minuscule λ ⟹ always simple | **proved** |
| checking the trivalent node suffices | **proved** (type E exhaustively) |
| type D closed form r + 2m ≤ 2, k ≤ 2 and spin | **proved** |
| type D closed form, k ≥ 3 | rests on Lemma S′ below |

The one open point is **Lemma S′**: in an extremal module, a special colour-t
line occurs only at a degree whose successor is already full. Verified on all
1990 colour-t graded pieces of D₄–D₇. Every graded submodule of I_κ that
violates upward closure has δ(v) = 2v_κ − vᵀCv > 0 — its dimension vector is not
extremal — so proving that implication would close the gap. See
`upclosure.py` and Remark 6.7 of the paper.

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
| `inwardlemma.py` | inward injectivity |
| `upclosure.py` | Lemma U, Lemma S′, the defect δ |

`RESULT.md`, `THEOREM.md`, `FINDINGS.md` are the working research logs, kept for
the record; they contain superseded formulations and are not the reference.
