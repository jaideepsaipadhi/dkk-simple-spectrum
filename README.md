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
r = |S ∩ {n−1, n}|. Then simple spectrum ⟺ **r + 2m ≤ 2**. This is now
unconditional; the route is described below.

**Reduction to one vertex.** Only the trivalent node t needs checking: if V_t is
multiplicity-free then so is V. Proved in both types, via inward injectivity on
I_κ — exhaustively in type E, and by a rank-uniform argument in type D.

## Quick start

    ./validate.sh           # the audit + all 27 cited scripts  (~5 min)
    ./validate.sh quick     # the audit alone                  (~2 min)
    ./validate.sh full      # ... and the whole exploration log (~45 min)

Exit status is 0 only if everything passes. Needs python3 (3.8+) and nothing
else — no packages, no network. Or run the pieces directly:

    python3 audit.py        # THE AUDIT: every link, plus negative controls
    python3 proof.py        # every step of the proof of the main theorem
    python3 verify.py       # DKK's printed posets, criterion vs solver, type D
    python3 sweepall.py     # aggregate criterion-vs-solver sweep (2179/2179)
    python3 inwardlemma.py  # the inward-injectivity lemma
    python3 master.py       # the staircase property, over all graded submodules
    python3 starchar.py     # (*o), the order-free form
    python3 reduce.py       # (*o) -> k explicit weights
    python3 canon.py        # tau, the canonical truncation, (T0)-(T2)
    python3 strands.py      # the strand decomposition (W1)-(W4)
    python3 thresh.py       # (T0)-(T2) and the arithmetic, D5-D60
    python3 localstar.py    # why (*o) cannot be localised further
    python3 rppmod.py       # the RPP hook product (Krylov's follow-up)
    python3 filltop.py      # the unconditional reductions (a) and (b)
    python3 claimG.py       # reduction (c), link by link
    python3 inwardproof.py  # inputs to the type-D proof of Lemma I
    python3 bigrade.py      # the bidegree is a function of (vertex, degree)

Each script prints a pass/fail line and exits 0/1. No dependencies beyond the
standard library. Linear algebra is over F_p; set `DKKP` to change the
characteristic (5, 7, 11, 13 all tested) and `DKKCAP` to raise the orbit cutoff
in `sweepall.py`.

## The manuscript

`paper.tex` / `paper.pdf` (20 pages).

## What is proved

The main theorem and the reduction to the trivalent node are unconditional and
uniform in type. The type-D closed form runs through the **staircase property**:

> If dim M is w-extremal then M is a staircase — top-truncated at the trivalent
> node t and an up-set at the two leaves.

That is what took the work. The chain that proves it, in order:

**1. Remove the word.** The property reduces to (\*): at every colour-t socle
step the partial module is an up-set at the chain neighbour t−1. Since the
partial module is itself extremal, (\*) follows from the order-free

> **(\*°)** an extremal M ⊆ I_κ with ⟨α_t^∨, λ − dim M⟩ ≥ 1 is an up-set at t−1

— 1998/1998 over D₅–D₁₀, and the side condition is sharp: every extremal module
failing to be an up-set at t−1 has ⟨α_t^∨, λ − dim M⟩ = 0, all 111 of them
(`starchar.py`). Extremality cannot be dropped: the same statement for arbitrary
submodules with (H) is **false**, first at κ = 4 (`localstar.py`).

**2. Remove the rank.** Let J = I ∖ {t−1, t}. A colour-j socle step for j ∈ J
grows the module only at vertex j, so "up-set at t−1" is constant on W_J-orbits,
and only the W_J-minimal weights survive — exactly the μ whose sole descent is
t−1, and there are exactly **k of them, for every rank**. 14373 instances
collapse to 68 representatives over D₅–D₁₀, k ≤ 6, the reduction step checked on
the modules themselves (1352/1352, `reduce.py`). At those weights V(μ_a) is the
greedy top truncation of I_κ at every vertex — statement **(B)**.

**3. Make it a closure check.** A top truncation is not determined by its
multiplicities: at a multiplicity-2 slot taking only 1, a *line* must be chosen.
τ, the leaf-exchanging automorphism, fixes I_κ and acts on such a slot with
eigenvalues ±1; the representatives always use the **+1 eigenline**. Define
T(v), the canonical top truncation, accordingly. Then dim T(v_a) = v_a, so by
rigidity **(B) ⟺ T(v_a) is a submodule** — no socle construction left
(`canon.py`: τ involutive 326/326, equivariant 552/552, T(v_a) = V(μ_a) 32/32).

**4. The strand decomposition (W).** At a chain vertex τ splits every graded
piece into eigenspaces, and the splitting is as clean as possible. With t = n−2,
so the top degree is 2t:

- **(W1)** each eigenspace is multiplicity-free — a 2-slot is one + line and one
  − line (8255/8255);
- **(W2)** a⁺_d(j) = 1 ⟺ |j−k| ≤ d ≤ j+k−2 and d ≡ j+k (mod 2), and
  a⁻_d(j) = a⁺_{2t−d}(j): the + strand is the Dirichlet half-line wave of the
  generator, the − strand its reflection in the top degree (27430/27430 against
  the profile, 8255/8255 against the eigen-dimensions);
- **(W3)** each strand is a contiguous run of degrees — + low, − high,
  overlapping in the multiplicity-2 window (2236/2236);
- **(W4)** every arrow restricts to an **isomorphism** between consecutive
  nonzero slots of the same strand (17706/17706).

(W2) follows from Lemma R: the complex there is τ-equivariant (at t the two leaf
summands are exchanged, so their sum is τ-stable), so in char ≠ 2 it splits, and
counting dimensions gives a recursion which, with a⁺_0(j) = δ_{jk} — the
generator e_κ is τ-fixed — determines everything. (W1) and (W3) are immediate
from (W2).

(W4) comes from the preprojective relation, which makes a dead arrow propagate:
with min(hi^ε_{j−1}, hi^ε_{j+1}) = hi^ε_j − 1, a death at (d,j) toward one chain
neighbour forces one at (d+1, other neighbour) toward j. The walk moves away
from the dead arrow: outward deaths reach vertex 1, which has a single chain
neighbour, contradicting Lemma R; inward deaths reach t and reduce to the
folded-leaf arrow being nonzero, which the relation at the leaf supplies.

**5. Everything about arrows becomes dimensions.** From (W):

- **(T0)** every arrow has rank min(dim source, dim target);
- **(T1)** a rank-1 arrow into a multiplicity-2 chain slot has image the
  τ-eigenline;
- **(T2)** the τ-eigenline is killed by exactly the arrows into 1-dimensional
  chain slots.

(T0) says the strand sets of adjacent slots are **nested**, and a crossing
collapses to an inequality violating j, k ≤ t or |i−j| = 1; (T1) and (T2) go the
same way. Checked from the closed form over D₅–D₁₀₀: 23,517,840 nesting
instances, 156,848 each of (T1), (T2).

**6. Pure arithmetic.** With (W), a τ-invariant submodule splits as N⁺ ⊕ N⁻ with
each part up-closed — and since (d,j) reaches (d+2,j), a subobject meets each
vertex in a *top interval*. So closure of T(v_a) is up-closedness of its two
strand-parts. Each strand has **min(j,k)** slots at vertex j, of which
**max(0, j+k−t−1)** are shared, so with u_j = v_j − min(j,k):

    s⁺_j = max(0, u_j, ⌈(u_j+o_j)/2⌉),   s⁻_j = v_j − s⁺_j

and up-closedness is exactly **0 ≤ s⁺_{j+1} − s⁺_j ≤ 1** and
**|s⁻_{j+1} − s⁻_j| ≤ 1**. Reading off the Bourbaki increments gives
u_{j+1} − u_j ∈ {−1,0,1}, negative only for j+1 ≤ a (both ends 0) or j+1 = t
(where o steps up, u+o is unchanged, and u_{t−1} ≤ o_{t−1}+1 stops a drop); and
where v jumps by 2 the plus-occupation rises by exactly 1, because either u
attains the maximum or o steps up. The leaf obligations close on three
identities at the node: v_t = 2k−2a−2, the truncation at t starts (and has its
threshold) at t−k+2a+2, and the leaves start at t−k+2a+3 — exactly one higher —
so every leaf obligation lands in a full slot.

**So the type-D closed form is unconditional.** End-to-end against the
independent quiver-variety solver: **5819/5819** over D₄–D₁₀, every κ.

## Machine verification

`audit.py` walks the whole chain in one run — 43 checks, ~90s — and includes
**negative controls**: statements that must *fail*, among them the one-neighbour
interval step that an earlier draft of the propagation lemma wrongly used. It
also re-runs the criterion at p = 5, 7, 11, 13, 1000003 (2,922,500 (T0)
instances and 164,160 (G1)/(G2) instances are checked from the closed form), and
verifies that the numbers quoted here and in the manuscript are the ones it
computes — the check that catches a stale count after a change.

    python3 audit.py          # the chain, negative controls, robustness
    DKKFULL=1 python3 audit.py   # ... and every script in the repo


`thresh.py` runs the whole arithmetic from closed forms over D₅–D₄₅, every k and
every a — the closed forms themselves checked against the actual modules first
(62/62 profiles, 9/9 dimension vectors, κ' = k):

    closure obligations, all of a resolved kind   6,225,317 / 6,225,317
    up-closed in each strand                      8,047,275 / 8,047,275
    occupation formula reproduces the truncation    459,282 / 459,282
    the two increment bounds                        445,178 / 445,178 each
    every leaf obligation lands in a full slot      446,982 / 446,982
    (P2) a drop only at p <= a or at t              148,994 / 148,994
    (G1),(G2) profile inequalities, to rank 80      164,160 / 164,160

`strands.py` checks (W1)–(W4) on the modules over D₅–D₁₇ and every κ, and
(T0)–(T2) from the closed form to rank 100. Wider single-purpose runs: the
strand counts to rank 120 (547,636), the threshold inequality to rank 100
(11,592,664), up-closedness to rank 60 (35,138,978).

## A negative result worth recording

The natural guess — that V is the greedy top truncation of I_κ at *every* vertex
— is **false**. It holds at 1185 of 1307 extremal modules over D₄–D₇, E₆, E₇, A₅;
the smallest failure is D₅, λ = ω₂, μ = (1,−1,0,0,0), where V₂ occupies degrees 2
and 6 while the top truncation would give 6 and 4. That element *does* have
simple spectrum, so the failure is not confined to degenerate cases. Top
truncation at t and at the two leaves does hold throughout.

Also: in D₄, w = s₂s₄s₃s₁s₂ is fully commutative but does **not** have simple
spectrum (v = (1,3,1,1); degree 2 carries colour 2 twice). With DKK's Example
5.12 this shows full commutativity fails in both directions.

## Krylov's follow-up: the RPP generating function

For every simple-spectrum w the RPP generating function of the poset is a hook
product. With a reduced word w = s_{i₁}⋯s_{i_ℓ}, inversion roots
β_p = s_{i₁}⋯s_{i_{p−1}}(α_{i_p}) and c-trace multiplicities c_p,

    RPP_{P(w)}(z) = ∏_p (1 − z^{ht β_p})^{−c_p}

— 342/342 over A₄–A₆, D₄–D₇, E₆, the hook multiset independent of the reduced
word (2181/2181). For λ-minuscule w every c_p = 1 and this is the
Peterson–Proctor product over Inv(w); beyond that locus (53 of the 342) the
product over Inv(w) has too few factors. Orientation matters: in the opposite
one it is a hook product only 110 times out of 342. See `rppmod.py`.
