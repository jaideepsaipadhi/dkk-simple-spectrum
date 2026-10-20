# Extremal torus fixed points and simple spectrum

Code and manuscript on the open problem of §5.4.2 of Dinkins–Karpov–Krylov,
*The quantum Hikita conjecture via quasimaps*,
[arXiv:2608.16746](https://arxiv.org/abs/2608.16746):

> "We consider it an interesting open problem to find a combinatorial
> characterization of such *w*."

**Scope.** DKK restrict neither the Dynkin type nor the weight (λ in their
Def. 5.1 is arbitrary). What is proved here answers the question **in type D,
for non-spin fundamental λ**. Type E at a non-minuscule fundamental weight, and
λ not fundamental in any type, are open and no criterion is claimed for them.
See §1.3 of the manuscript, which states the boundary precisely.

## The manuscript

`typeD.tex` / `typeD.pdf` (44 pages). This supersedes `paper.tex`, which is
retained only as a record of an earlier and partly superseded approach — see
"Superseded files" below before reading it.

## The result

Let Γ be a simply-laced Dynkin diagram, λ = ω_κ a fundamental weight, w ∈ W of
minimal length in wW_λ, and μ = wλ. Then M(λ−μ, λ) is a point with a unique
T-fixed representation V. Let Π be the preprojective algebra of Γ and I_κ the
injective envelope of the simple S_κ.

**Main theorem** (uniform in type, characteristic-free). Read a reduced word for
w from the right, recording at each step c = ⟨α_i^∨, μ_current⟩; let σ be the
resulting sequence of vertices with multiplicities. Then

    V  ≅  soc_σ(I_κ)

as graded modules, where soc_σ is the iterated socle. Simple spectrum ⟺ no two
elements of the resulting coloured graded set share a colour and a degree.

**Type D closed form.** For λ = ω_k with 2 ≤ k ≤ n−2, over an arbitrary field,

    simple spectrum  ⟺  v_t ≤ 2  ⟺  r + 2m ≤ 2

where t = n−2 is the trivalent node, v = λ−μ, and writing μ = Σ_{i∈S} η_i e_i in
Bourbaki coordinates, m = #{i ≤ n−2 : i ∈ S, η_i = −1} and r = |S ∩ {n−1, n}|.
The two forms agree because v_t = r + 2m. For the three spin weights
ω_1, ω_{n−1}, ω_n every such w is λ-minuscule and the spectrum is always simple;
that is a genuinely different statement, since v_t is unbounded there.

**A necessary condition at the node.** If the c-trace exceeds 1 at any
occurrence of the trivalent node, the spectrum is not simple. Proved for every
fundamental λ in type D. The same conclusion holds in E₆ as a computation: of
the 1272 extremal weights of the six fundamental λ, the 831 whose c-trace
exceeds 1 at the trivalent node all fail to have simple spectrum. No proof of
the E₆ statement is known to us.

## Quick start

    ./validate.sh           # the audit + all cited scripts
    ./validate.sh quick     # the audit alone

Exit status is 0 only if everything passes. Needs python3 (3.8+) and nothing
else — no packages, no network. The scripts for the current results:

    python3 socle.py          # the extremal module as an iterated socle
    python3 socval.py         # soc_sigma cross-validated against a direct solver
    python3 verify_nodeone.py # the node condition, D4-D6 exhaustive
    python3 verify_nodeone.py 8 --fast
    python3 efeat.py E6       # the type-E criterion search (1272 weights, 100s)
    python3 tricheck.py       # the node condition across type D

Linear algebra is over F_p; set `DKKP` to change the characteristic.

## Ranges of the checks

Every statement of the manuscript is proved there or cited; nothing in its proof
rests on computation. The computations are independent confirmations, tabulated
in §14.2 of `typeD.pdf`. The largest:

| statement | range | instances |
|---|---|---|
| the criterion vs a direct solver | full orbits D₄–D₁₁, p = 2,3,5 | 241008 per p |
| specialisation to p = 2, slotwise | full orbits D₄–D₁₀, all k | 77196 weights |
| sharpness of the side condition in (\*°) | full orbits D₄–D₁₀, all k | 77196 modules, 2778 non-up-sets, all with c_t = 0 |
| the ℤ-form of a chain slot | D₄–D₂₂, all k | 7315 multiplicity-2 slots |
| the block induction δ_p = c(p) | full orbits D₄–D₈ | 519184 blocks |

The criterion check was produced by a from-scratch implementation sharing no
code with the rest of the project.

## Negative results worth recording

**Full commutativity is neither necessary nor sufficient.** In D₄,
w = s₂s₄s₃s₁s₂ is fully commutative and does *not* have simple spectrum
(v = (1,3,1,1); degree 2 carries colour 2 twice). DKK's Example 5.12 gives a w
with simple spectrum that is not fully commutative.

**The greedy top truncation at every vertex is false.** The smallest failure is
D₅, λ = ω₂, μ = (1,−1,0,0,0), where V₂ occupies degrees 2 and 6 while the top
truncation would give 6 and 4. That element *does* have simple spectrum, so the
failure is not confined to degenerate cases. Top truncation at t and up-sets at
the two leaves do hold throughout.

**The larger parabolic does not work.** With J = I ∖ {t−1, t} the side condition
of (\*°) is not W_J-invariant, so the reduction through J — which would give k
representatives — is not closed. The reduction runs through
J' = {1, …, t−2} instead and produces **9(k−1)** representatives. An earlier
draft of this repository claimed k; that was wrong, and §14.3 of the manuscript
gives the counts behind the correction.

**Simple spectrum is not characterised by the c-trace shape.** The guess that it
is equivalent to a c-trace of the form (1,…,1,3) is false: of the nine elements
of E₆ ∪ E₇ with a c-trace entry 3 and simple spectrum, five have other shapes.

## Krylov's follow-up: the RPP generating function

For the simple-spectrum w tested, the RPP generating function of the poset is a
hook product: with inversion roots β_p = s_{i₁}⋯s_{i_{p−1}}(α_{i_p}) and
c-trace multiplicities c_p,

    RPP_{P(w)}(z) = ∏_p (1 − z^{ht β_p})^{−c_p}

over A₄–A₆, D₄–D₇, E₆, with the hook multiset independent of the reduced word.
For λ-minuscule w every c_p = 1 and this is the Peterson–Proctor hook-product
formula over Inv(w); beyond that locus the product over Inv(w) has too few
factors. See `rppmod.py`.

**Prior art.** This is not new territory. Peterson–Proctor proves every
d-complete poset is a hook length poset; Proctor classified the d-complete
posets; and "hook length poset" is the established term, with the relevant
structure theory in Ishikawa–Tagawa, *Leaf posets and multivariate hook length
property*. The skew case is Naruse–Okada's skew hook formula, quoted by DKK
themselves at the end of their §5.3. Anyone extending this should read those
first.

## Superseded files

These are kept for the record and are **not** current. Where they disagree with
`typeD.pdf`, the manuscript is right.

- `paper.tex` / `paper.pdf` — earlier manuscript, superseded by `typeD.tex`.
- `RESULT.md`, `THEOREM.md` — earlier write-ups. Both predate the correction of
  the representative count from k to 9(k−1), and `RESULT.md` claims the open
  problem of §5.4.2 without the type-D / fundamental-λ qualification stated
  above.
- `FINDINGS.md` — working notes from the computational attack.

In particular, the commit message `"Problem 1 solved unconditionally"` in this
repository's history refers to the ground field — the type-D criterion holds in
every characteristic, p = 2 included — and **not** to the problem of §5.4.2 in
general. See the Scope note at the top.
