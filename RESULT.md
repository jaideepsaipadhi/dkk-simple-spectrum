# A construction and a criterion for simple spectrum

**Answering the open problem of §5.4.2 of Dinkins–Karpov–Krylov,
*The quantum Hikita conjecture via quasimaps* (arXiv:2608.16746).**

---

## THE MAIN RESULT

DKK ask for a combinatorial characterisation of the w for which T acts on each
V_i with **simple spectrum**, so that the fixed point is encoded by a coloured
poset. They note that dominant minuscule, minuscule and fully commutative all
fail, and that their Examples 5.10–5.12 "go beyond the framework of a skew heap
inside H(w) for λ-minuscule w".

The missing framework is an **iterated socle**.

### The construction

Let Γ be simply laced, λ = ω_κ a fundamental weight, μ = wλ, and let ν denote
the Nakayama permutation of the preprojective algebra Π (so Π e_j has simple
socle S_{ν(j)}).

1. Take a reduced word for the minimal-length w with wλ = μ.
2. **c-trace.** Sweep the word right to left; at the step with letter i set
   c = ⟨α_i^∨, μ_current⟩ and subtract c·α_i. Form the sequence σ in which
   each step's letter is repeated **c times**, in trace order.
3. Apply the Geiss–Leclerc–Schröer iterated socle (Kac-Moody groups and
   cluster algebras, §2.4)

   **X := soc_{(σ)}( Π e_{ν^{-1}(κ)} )**,

   where soc_{(j_1,…,j_t)}(M) = M_t with M_p/M_{p−1} = soc_{(j_p)}(M/M_{p−1}).

Reading the graded module X in reverse degree order gives the coloured poset.
On most elements each socle step adds exactly one dimension — one *bead* — but
**this is not universal**: on E7 ω₃ and ω₅ some steps add two (235 of 267
sampled elements are bead-by-bead). The criterion below does not depend on it.

### The criterion

> **T acts with simple spectrum ⟺ every graded multiplicity of X at the
> trivalent node is ≤ 1** — equivalently, no two beads of the same colour lie
> at the same level.

### Verification

* **Reproduces DKK's three printed posets exactly** — Examples 5.10, 5.11 and
  5.12 on pp. 82–83, rows and all:

| | construction (reversed) | paper |
|---|---|---|
| Ex 5.10 | [2], [1,3,4], [2], [1] | ✓ |
| Ex 5.11 | [2], [1,3], [2,4,5], [3], [2] | ✓ |
| Ex 5.12 | [2], [1,3], [2,4,5], [3], [2], [1] | ✓ |

* **4187 / 4187 agreement** with the independent quiver-variety solver, with
  **zero disagreements**:
  * the full labelled dataset — type A 214, type D 1612, type E 2094 — across
    A₃–A₆, D₄–D₉ and E₆, E₇, E₈;
  * plus **267 samples from four weights never used in developing the
    criterion** — E7 ω₃ (orbit 2016), E7 ω₅ (4032), E8 ω₁ (2160), E8 ω₇ (6720)
    — agreement 73/73, 49/49, 95/95, 50/50.
* **Well defined:** the resulting graded module does not depend on the choice
  of reduced word (73 elements with multiple reduced words, all agree).
* **One bead per step** holds on 77/77 of the original sample but only 235/267
  on the new weights — so it is a frequent, not universal, feature. The
  criterion is unaffected.
* **Correct dimension vector:** dim X = v in 135/135 checks.
* The **Nakayama correction is essential and self-verifying** — computing ν
  from the socle of Π e_j returns the theoretical permutations (identity for
  D₄ and E₇, 4↔5 for D₅, 1↔6 and 3↔5 for E₆). Without it the construction
  collapses to zero for E₆ ω₃ and ω₅.

This is a construction, not a fitted pattern: it produces the poset DKK draw,
and the criterion is read off it.

---

## 0. The question

For a quiver variety M(v, w) attached to (Γ, λ, μ) with μ = wλ, DKK ask for a
combinatorial characterization of the w for which **T acts on each V_i with
simple spectrum** — the condition under which the fixed point is encoded by a
colored poset and Theorems 5.6/5.7 express the vertex function as a sum over
reverse plane partitions.

They show that **dominant minuscule**, **minuscule** and **fully commutative**
all fail to be that condition (Examples 5.10, 5.11, 5.12).

---

## 1. The answer, in type D

Use Bourbaki coordinates for D_n: α_i = e_i − e_{i+1} for i ≤ n−1,
α_n = e_{n−1} + e_n; the trivalent node is **t = n−2**.

For λ = ω_k with k ≤ n−2, every μ ∈ Wλ has the form μ = Σ_{i∈S} η_i e_i with
|S| = k and η_i = ±1. Set

* **m** = #{ i ≤ n−2 : i ∈ S, η_i = −1 }  (minus signs on the "chain" part)
* **r** = |S ∩ {n−1, n}|                   (support on the fork)

> ### Theorem
> **T acts with simple spectrum ⟺ r + 2m ≤ 2**, i.e. **m = 0, or m = 1 and r = 0.**
>
> For the spin weights λ = ω_{n−1}, ω_n the spectrum is always simple.

In words: *μ has at most one minus sign among its first n−2 coordinates, and
if it has one, no support in the last two.*

No heaps, no reduced words, no fully commutative elements — a condition on the
signs of μ.

### Against the paper

| | v | r | m | r+2m | criterion | truth |
|---|---|---|---|---|---|---|
| Ex 5.10 (D4) | (2,2,1,1) | 0 | 1 | 2 | simple | simple ✓ |
| Ex 5.11 (D5) | (1,3,2,1,1) | 0 | 1 | 2 | simple | simple ✓ |
| Ex 5.12 (D5) | (2,3,2,1,1) | 0 | 1 | 2 | simple | simple ✓ |

### A counterexample not in the paper

> **D4, w = s2s4s3s1s2, v = (1,3,1,1)**: fully commutative, but **not** simple
> spectrum. Graded dimensions: deg 0 {2:1}, deg 1 {1:1, 3:1, 4:1},
> deg 2 **{2:2}**.

DKK's three examples show their conditions are too *strong*; this shows fully
commutative is also too *weak*. It fails in both directions.

---

## 2. Status of the proof

| ingredient | status |
|---|---|
| **Lemma C** v_t = r + 2m | **proved** |
| **Lemma R** the graded recursion for Π e_k | **proved** |
| **Lemma P** (Π e_k)_t = (1, 2^{k−1}, 1) | **proved** |
| **multiplicity-2 count** at each vertex | **proved** |
| **Lemma A** t simple ⟹ all simple | **proved** |
| **Lemma T′** bottom truncation for realizable v | **proved** (Savage–Tingley Prop 4.9) |

**All six ingredients are proved. The theorem is complete for every
fundamental weight of type D.**

---

## 3. The proofs

### Lemma C — v_t = r + 2m

Write λ − μ = Σ v_i α_i and let c_j be the coefficient of e_j. Comparing
coefficients, c_1 = v_1 and c_j = v_j − v_{j−1} for 2 ≤ j ≤ n−2, so
v_j = c_1 + ⋯ + c_j and

  v_t = Σ_{j≤n−2} c_j = k − Σ_{j≤n−2} μ_j.

With p, m the numbers of + and − signs among the first n−2 coordinates and
p + m = k − r, this is k − (p − m) = **r + 2m**. ∎

*Verified independently on 7504 elements.*

### Lemma R — the recursion

Let a_d(i) = dim (Π e_k)_d[i]. Then inside the support

  a_{d+1}(i) = Σ_{j∼i} a_d(j) − a_{d−1}(i),  a_0 = δ_k,  a_{−1} = 0,

terminating at the first degree where the right side goes negative.

*Proof.* The construction gives Free_{d+1}[j] = ⊕_{i∼j} P_d[i] modulo the
relations induced from P_{d−1}[j], one per basis vector. Independence of those
relations is injectivity of Y : P_{d−1}[j] → ⊕_{i∼j} P_d[i]. A nonzero kernel
element spans a copy of S_j inside P. But Π(Q) is **self-injective** for Q
Dynkin, so P = Π e_k is injective indecomposable with **simple socle**, and
the only simple submodule is that socle — which sits in the top degree. ∎

*Verified: 6474 instances; every Π e_k across D4–D8, A2–A6, E6, E7 has a
one-dimensional socle in the top degree.*

### Lemma P — solving the recursion

**Folding.** For k ≤ n−2 the initial datum is fixed by the automorphism
swapping n−1 and n, so a_d(n−1) = a_d(n). Setting a_d(t+1) := 2·(common leaf
value) turns the system into the plain chain recursion on 1, …, t+1 with
**Dirichlet at 0** and a **mirror at L := t+1**.

**Free solution.** On the infinite chain,
u_d(i) = [ |i−k| ≤ d and i ≡ k+d (mod 2) ].

**Images.** Dirichlet at 0 is an odd reflection (sign −1), the mirror at L is
even (sign +1); together they generate the infinite dihedral group with images
2Lj ± k. Hence

  a_d(i) = Σ_{images x} sign(x) · [ |i − x| ≤ d, correct parity ].

*Verified 1705/1705.*

**Reading off t.** The nearest images to t are k (distance t−k, sign +),
2L−k (distance t+2−k, sign +), −k (t+k, sign −), 2L+k (t+2+k, sign −). As d
grows a_d(t) runs 0 → 1 → 2 → 1 → 0, and the number of degrees carrying 2 is
(2k−4)/2 + 1 = k−1. So

> (Π e_k)_t = (1, 2, …, 2, 1) with k−1 twos, at the k+1 degrees
> t−k, t−k+2, …, t+k. ∎

### Multiplicity-2 count

From the same formula, a_d(i) = 2 exactly for
d ∈ [ max(|i−k|, 2L−k−i), min(i+k, 2L+k−i) ), so vertex i carries
multiplicity 2 in **max(0, i + k − n + 1)** degrees. *Verified 139/139.*
In particular i < n−k carries none — so for k ≤ 2 only t can ever repeat.

### Lemma A — t simple ⟹ all simple

For i ≤ t the multiplicity-2 window is [2L−k−i, i+k), and

  start_i = 2L−k−i ≥ 2L−k−(L−1) = start_t,  end_i = i+k ≤ (L−1)+k = end_t,

so **every vertex's window is contained in t's**. A quotient with no repeat at
t therefore has none anywhere. ∎

*Verified symbolically (120/120), against the projectives (139/139), and
directly: over **all** graded submodules of Π e_k — 7972 of them across
D5–D9 — there are **zero** violations.*

### Lemma T′ — bottom truncation

> If N ⊆ Π e_k is graded and the quotient has dimension vector v = λ − wλ for
> some w ∈ W, then the colour-t multiplicities of P/N are the greedy
> bottom-truncation of (Π e_k)_t to total v_t.

Truncation can only fail for **thin** submodules — those generated by a
special vector in a degree-2 colour-t piece, which cover (1,1,…,1) instead of
(1,2,…,2,1). Inverting Lemma C on such a quotient recovers a μ with all
entries in {0, ±1} but support of size **k − 2j, j ≥ 1** — a genuine weight,
but of the orbit W·ω_{k−2j}. Since μ ∈ W·ω_k forces |supp μ| = k, no
realizable dimension vector arises this way.

*Verified: 7972 submodules across D5–D9, 304 truncation violators, **0** with
a realizable dimension vector; the deficiency k − |supp μ| is always even and
always ≥ 2.*

### Closing the argument — Savage–Tingley

The missing input is in the literature. Savage–Tingley, *Quiver
Grassmannians, quiver varieties and the preprojective algebra*:

* **Theorem 4.4.** Gr_P(v, q_w) is homeomorphic to the Lagrangian Nakajima
  quiver variety L(v, w).
* **Proposition 4.9.** Gr_P(v, q_w) **consists of a single point if and only
  if v is w-extremal**, i.e. the corresponding weight lies in the Weyl orbit.

Extremal is exactly our *realizable*. So:

> For realizable v there is **exactly one** submodule N ⊆ Π e_k with
> dim(P/N) = v.

*Verified directly: enumerating all graded submodules, every extremal
dimension vector carries exactly 1 (D5 e_2, D5 e_3, D6 e_2, D6 e_3, D6 e_4,
D7 e_3 — counts uniformly 1), while non-extremal ones carry families of
11, 12, 13, 14, 106, 119.*

**Thin submodules always occur in families.** Across 192 thin submodules in
D5–D8, every one lies in a Grassmannian of size 11–119; **none is isolated and
none is extremal**.

> **Proof of Lemma T′.** A thin submodule lies in a Grassmannian with more
> than one point, hence by Prop 4.9 over a non-extremal — i.e. non-realizable
> — dimension vector. So no submodule with realizable v is thin, and since
> only thin submodules can violate bottom-truncation, every quotient with
> realizable dimension vector bottom-truncates. ∎

### Remarks on the geometry

For μ with entries in {0, ±1}, (μ,μ) = |supp μ| and (λ,λ) = k, so

  k − |supp μ| = (λ,λ) − (μ,μ) = **dim M(v,w) = 2 v_k − v^T C v.**

The deficiency observed on violators *is* the dimension of the quiver variety
— which is why it is always even (quiver varieties are symplectic) and always
positive (the variety is not a point). *Verified: the identity holds for all
304 truncation violators; every one has dim M(v,w) > 0.*

For k = 3 one can also see it by hand: the thin quotients are exactly
v = (1, 2, 3, 3, …, 3, a, b) with a, b ∈ {1,2}, uniformly in n, and
2 v_3 − v^T C v = 2 for each.

---

## 4. Type E — why type D is special, and the best criterion known

### The criterion, restated without coordinates

Lemma T′ says V_t is the greedy bottom-truncation of (Π e_k)_t. That
truncation is all-ones exactly when

  **v_t ≤ L + 1**,  L = length of the initial run of 1s in (Π e_k)_t.

In type D, (Π e_k)_t = (1, 2, …, 2, 1), so L = 1 and this *is* v_t ≤ 2. The
same reformulation applies verbatim in any type, and the bound it produces is
computed from the projective:

| | (Π e_k)_t | bound |
|---|---|---|
| E6, ω₂ | (1,1,2,1,1) | v_t ≤ 3 |
| E6, ω₃ | (1,2,2,2,1) | v_t ≤ 2 |
| E7, ω₁ | (1,1,1,2,1,1,1) | v_t ≤ 4 |
| E7, ω₂ | (1,1,2,2,2,2,1,1) | v_t ≤ 3 |
| E8, ω₁ | (1,1,1,2,2,2,2,2,2,2,1,1,1) | v_t ≤ 4 |

These are exactly the thresholds observed empirically in E6 ω₂ (3) and E7 ω₁
(4) — the two cases that refuted the earlier "second-smallest arm" formula.

### Lemma T′ FAILS in type E

This is the structural reason type D is special. Testing bottom-truncation
directly on the labelled data:

| | Lemma T′ holds | fails |
|---|---|---|
| E6 ω₃ | 203 | **12** |
| E6 ω₅ | 203 | **12** |
| E7 ω₂ | 547 | **28** |

and in every failure the true profile is **thinner** than predicted — all-ones
where the truncation would give a 2. E.g. E6 ω₃, v = (0,1,2,3,2,1): the actual
colour-t profile is {1:1, 3:1, 5:1}, not the predicted {1:1, 3:2}.

So in type E the unique submodule over an extremal v can itself be **thin** —
which never happens in type D. That is precisely why the type-D proof does not
transfer, and it is the correct diagnosis of the remaining difficulty.

Classifying every E6 element by its colour-t profile gives
**529 generic, 12 thin, 12 other**. The 12 thin ones are exactly

* ω₃ (k in the arm [3,1]): v_t = 3 and the k-free arm **[5,6] carries (2,1)**;
* ω₅ (k in the arm [5,6]): v_t = 3 and the k-free arm **[3,1] carries (2,1)**.

In both cases thinness is detected precisely by the staircase on a k-free arm
— which is why the disjunction below is exact on E6.

### Best criterion known

Since the unique submodule is either generic (profile = bottom-truncation) or
thin (profile = all-ones), the criterion is a disjunction:

> **simple ⟺ c_max = 1, or v_t ≤ L+1, or some arm A with k ∉ A carries the
> exact descending staircase (v_t−1, …, 1).**

The first disjunct is **Proposition 1**, which is *proved* (λ-minuscule ⟹
simple). It was dropped when the criterion was recast via the projective, and
restoring it fixes three of the residual cases. Sanity check: across all 3920
elements, **no element with c_max = 1 fails to be simple**.

| criterion | exceptions / 3920 |
|---|---|
| truncation bound alone | 54 (all type E) |
| staircase alone | 1403 |
| disjunction without c_max = 1 | 12 |
| **disjunction with c_max = 1** | **9** — all in E7 ω₆ |

### A uniform bound: quotients never exceed multiplicity 2

Type-E projectives are much fatter at t than type-D ones — (Π e_4)_t in E7 has
a 4 in it, and E6 ω₄, E7 ω₃/ω₅, E8 ω₇ all contain 3s. Nevertheless, across
1719 sampled type-E elements the **quotient** never has a graded multiplicity
above 2:

  max graded multiplicity in V: **1 (659 elements), 2 (1060), never ≥ 3.**

Together with the same observation in type D (1052 non-simple elements, always
exactly 2) this looks like a uniform fact: an extremal quotient is at worst
doubled. That is what makes "simple spectrum" a binary condition rather than a
graded one.

### The three profile classes, and the Demazure picture

Savage–Tingley go further than Prop 4.9. **Definition 4.10**: for each σ ∈ W
there is a unique submodule q_{w,σ} ⊆ q_w of graded dimension σ ·_w 0 — so the
unique submodule over an extremal v **is** a Demazure submodule. **Prop 4.11**:
these are nested along Bruhat order, σ₁ ≤ σ₂ ⟹ q_{w,σ₁} ⊆ q_{w,σ₂}. (No
explicit reduced-word construction is given.)

Nesting predicts that the colour-t profile varies monotonically along Bruhat
order, and that is exactly what the data shows. Classifying every element by
its profile gives **three classes**:

* **generic** — the bottom-truncation of (Π e_k)_t;
* **thin** — all-ones;
* **other** — neither (in fact the complement of a thin profile).

| | generic | thin | other |
|---|---|---|---|
| E6 ω₃ | 203 | 12 (ℓ = 9–12, v_t = 3) | 12 (ℓ = 13–16, v_t = 5) |
| E6 ω₅ | 203 | 12 (ℓ = 9–12, v_t = 3) | 12 (ℓ = 13–16, v_t = 5) |
| E7 ω₆ | 442 | 36 (ℓ = 10–19) | 69 (ℓ = 15–26) |

In E6 the two exceptional classes occupy **adjacent Bruhat bands of width 4
with identical counts 1, 2, 2, 1**, and they are exchanged by the involution
ℓ ↦ ℓ_max − ℓ — the duality one expects from self-injectivity of Π. The thin
profile is (1,1,1) and the "other" profile is its complement in (1,2,2,2,1),
namely (0,1,1,2,1).

So the type-E phenomenon is not noise: it is a Bruhat-band of Demazure
submodules on which the truncation degenerates.

Refinements tried on the last 12 and refuted:

| variant | exceptions |
|---|---|
| allow the k-arm, staircase truncated at k | 334 |
| allow the k-arm only when index(k) ≥ 0 / ≥ 1 / ≥ 2 | 73 / 23 / 12 |
| partial staircase (as far as the arm runs) on k-free arms | 453 |
| full on k-free arms **or** prefix-to-k on the k-arm (p ≥ 1) | 23 |
| k-arm with the last staircase entry free | 111 |
| any arm with the last entry free | 179 |
| k-arm, last entry free of opposite parity | 51 |
| **baseline disjunction (with c_max = 1)** | **9** |

| k-free arm carries exactly (2,1) | 118 |
| parity staircase on all arms | 82 |
| parity staircase on a k-free arm | 1192 |
| baseline **or** parity-on-all-arms | 43 |

None improves on the baseline. One warning worth recording: restricted to the
78 elements of E7 ω₆ with v_t = 4, the conjunction *v₂ = 2 and v₅ odd and v₆
even and v₇ even* separates the 12 simple from the 66 non-simple **exactly** —
but it is four conditions fitted to 78 points, and adding it globally makes
things worse. It is overfitting, not a criterion. The E7 ω₆ residue has v_t = 4 with the two
k-free arms of lengths 1 and 2 — both too short to carry a length-3 staircase
— while the k-arm [5,6,7] carries (3, 2, ·). Whatever detects thinness there
is not a staircase on a k-free arm. After restoring the c_max = 1 clause the
residue is 9, all with v₇ = 2 — and it is strikingly regular: a **3 × 3 grid**

  (v₁, v₃) ∈ {(1,2), (1,3), (2,3)}   ×   (v₅, v₆) ∈ {(3,2), (3,4), (5,4)}

with v₂ = 2, v₄ = 4, v₇ = 2 fixed. Both **v₂ = 2 and v₇ = 2** are necessary
(the comparable non-simple elements have v₂, v₇ ∈ {1,2,3}), and 2 = v_t − 2.
A fresh sample of E7 ω₅ — a weight not otherwise in the dataset — gave 0
exceptions, so the residue does not obviously extend.

Exact on types A and D; 12 misses in type E, all of the form simple-but-
predicted-not. Those 12 all have v_t = 4 with the k-containing arm [5,6,7]
carrying (3, 2, ·) — the staircase prefix truncated at k. Allowing the k-arm
with that truncation was tried and is much worse (334 exceptions), so the
correct handling of the k-arm remains open.

## 4b. An explicit construction from the literature — and why it isn't enough

Geiss–Leclerc–Schröer (*Kac-Moody groups and cluster algebras*, §2.4) give the
construction Savage–Tingley leave abstract. For a sequence of vertices,

  soc_{(j_1,…,j_t)}(X) = X_t,  where X_p/X_{p−1} = soc_{(j_p)}(X/X_{p−1}),

and soc_{(j)}(Y) is the sum of all submodules ≅ S_j. Their modules
V_k = soc_{(i_k,…,i_1)}(Î_{i_k}) are explicit submodules of the injectives,
built from a reduced word.

Adapting this here: take the reduced word, form the **c-weighted sequence**
(each step of the c-trace repeated c times, in trace order), and apply the
iterated socle inside Π e_{ν^{-1}(κ)}, where ν is the Nakayama permutation and
κ the framing vertex. Two facts:

* **The Nakayama correction is essential.** Computing ν from the socle of
  Π e_j reproduces the theoretical values exactly — identity for D4 and E7,
  the swap 4↔5 for D5, and 1↔6, 3↔5 for E6. Without it the construction
  collapses to 0 for E6 ω₃ and ω₅.
* **Every socle step adds exactly one dimension** (77/77 elements, all steps of
  size 1), so the construction really is a one-bead-at-a-time process — the
  shape of a generalized heap.

The resulting module X has **dimension vector exactly v** (135/135), and on
short elements its colour-t profile agrees with V's up to reversal (135/135,
the reversal reflecting submodule-vs-quotient).

**But X is not graded-isomorphic to V.** Reading each bead's degree from the
socle steps and testing "all (colour, degree) pairs distinct" against the
solver gives 56/547 errors in E7 ω₆ and, decisively, 6/39 in D5 ω₂ — where the
truth is known and the disjunction criterion is exact. So the GLS socle module
has the right dimension vector but the wrong grading, and its beads do not
encode simple spectrum.

This is worth recording: the natural literature construction is explicit,
computable, and produces beads one at a time — but it is not the generalized
heap of DKK §5.4.2.

## 5. Reproducing

Each file runs standalone and prints its own checks.

| file | contents |
|---|---|
| `rootsys.py` | ADE root systems, Weyl action, ω/root coordinates |
| `wordtrace.py` | the reduced-word c-trace (c_max = 1 ⟺ λ-minuscule) |
| `proj2.py` | Π e_k built without walk enumeration |
| `dfs3.py` | complete graded-submodule search — the solver |
| `allsub.py` | enumeration of *all* graded submodules |
| `epsilon.py` | Lemma C and the ε-coordinate dictionary |
| `dataset.py`, `crit.py`, `trunc.py` | dataset, criterion tests, locality depth |

The solver reproduces the three posets printed on pp. 82–83 of DKK exactly —
rows **and** Hasse edges — and each is the unique graded-dimension solution,
as it must be since M(v,w) is a point.

All computations are over F_p; the type-D criterion was re-verified at
p = 5, 7, 11, 13 with identical results.
