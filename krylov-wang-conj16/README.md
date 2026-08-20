# Searching for counterexamples to Conjecture 1.6 (Krylov–Wang)

Computational attack on Conjecture 1.6 of Krylov–Wang, *A formula for the
q-character of functions on the nilpotent cone of some Lie algebra
representations*, [arXiv:2608.03314](https://arxiv.org/abs/2608.03314).

> **Conjecture 1.6.** Let $G$ be semisimple. If $V$ is a Hesselink-type
> representation, then $V$ is cofree. Moreover, (2) and Theorem 1.2 hold for $V$.

## What this computes

**The Hesselink-type test (exact).** $V$ is Hesselink-type when

    X_V = [ sum_w (-1)^l(w) t^{w rho} prod_{alpha in S+} (1 - q t^{-w alpha}) ]
          / [ sum_w (-1)^l(w) t^{w rho} ]

is a `Z[q]`-multiple of the trivial character. Expanding the product over
subsets `T` of `S+` and writing `sigma_T` for the sum of the weights in `T`,
the numerator is `sum_T (-q)^{|T|} A(rho - sigma_T)` with
`A(mu) = sum_w (-1)^l(w) t^{w mu}`. Each `A(rho - sigma_T)` either vanishes
(when `rho - sigma_T` is W-singular) or equals `+-A(nu)` for a dominant regular
`nu`, so by the Weyl character formula

    X_V = sum_T (-q)^{|T|} * sign(T) * chi_{nu(T) - rho}.

This is a **finite, exact** computation — no series truncation — and the
subsets are summed by convolution over the weights of `S+`, so the cost is
`|S+|` times the number of distinct partial sums rather than `2^{|S+|}`.
`V` is Hesselink-type iff every coefficient except the one at the trivial
character vanishes.

**The discriminator: identity (2).** For each Hesselink-type `V` the script
independently computes both sides of

    X_V = dim_q C[V^{rho^vee}]^h / dim_q C[V]^G

as power series. `dim_q C[V]^G` is obtained by multiplying
`ch_q C[V] = prod_{alpha in S} 1/(1 - q t^{-alpha})` by the Weyl denominator
and reading the coefficient at `t^rho`; `dim_q C[V^{rho^vee}]^h` counts
functions `m : S0 -> Z_{>=0}` with `sum m(a) a = 0`, graded by `|m|`.

**A Hesselink-type `V` that fails identity (2) refutes Conjecture 1.6.**

**The stronger test: Theorem 1.2 directly.** Proposition 2.1 derives Theorem
1.2 from (2) *only for cofree `V`* — it goes through Lemma 2.3,
`C[N_V] (x) C[V]^G = C[V]`, which needs cofreeness. So a Hesselink-type `V`
that is **not** cofree can satisfy (2) and still violate Theorem 1.2, and the
identity-(2) test would never see it. `nullcone.py` therefore builds

    C[N_V] = C[V] / (C[V]^G_+ C[V])

directly by linear algebra — invariants as the weight-zero polynomials killed
by every raising operator, then the ideal they generate, degree by degree —
decomposes each graded piece into irreducible characters, and `mq.py` compares
it against `sum_lambda M_q(lambda) chi_lambda`. This is implemented for
`G = SL2^k`, where the representations are explicit and where §3.4 places the
classification question. Run it with `sweep12.py`.

The detector is calibrated in both directions before every sweep: `V_1`,
`V_1+V_1`, `V_2` must agree (Theorem 1.2 is proved there) and `V_3`, `V_4`
must disagree (§2.1 gives them a modified `P`). It correctly reports exactly
that.

## Running it

No dependencies beyond the Python standard library (`python3` ≥ 3.8).

```bash
python3 hess.py       # reproduces the paper's Examples 4.2 and 4.3, and the adjoint cases
python3 series.py     # validates identity (2) on the cases where Thm 1.2 is proved
python3 mq.py         # validates the direct Theorem 1.2 test on known cases
python3 census.py     # search via identity (2)   -- broad, cheap
python3 sweep12.py    # search via Theorem 1.2    -- narrower, decisive
```

`sweep12.py` options:

```bash
python3 sweep12.py --k 3 --maxdim 12 --maxdeg 5
python3 sweep12.py --k 4 --maxdim 16 --maxdeg 4 --out sl2fourth.json
```

Cost is dominated by the nullcone computation and grows quickly in both
`dim V` and `--maxdeg`: degree `d` needs `binomial(dim V + d - 1, d)`
monomials. `dim V = 10` at `--maxdeg 4` is a few seconds per representation.
Results stream to stdout and are written incrementally, so a run can be
interrupted without losing work.

`census.py` options:

```bash
python3 census.py --groups A1 A1^2 A2 A1^3 B2 G2 A1xA2 A3 \
                  --maxdim 12 --maxparts 3 --order 8 --out census.json
```

- `--groups` — group names: `A1`, `A2`, …, `B2`, `B3`, `C3`, `D4`, `G2`, `F4`,
  and products written with `^` or `x`: `A1^3`, `A1xA2`, `A1^2xB2`.
- `--maxdim` — bound on `dim V`.
- `--maxparts` — maximum number of irreducible summands of `V`.
- `--order` — truncation order in `q` for the identity-(2) test.

Representations on which some simple factor acts trivially are skipped: they
are really representations of a smaller group and are covered by that group's
own sweep. Adding trivial summands changes neither `C[N_V]`, `X_V`, nor
cofreeness, so nothing is lost.

**Where to push.** `X_V` is cheap; the identity-(2) test is the expensive part
and only runs on Hesselink-type representations, which are rare. So raising
`--maxdim` is usually affordable, and that is the direction with the most room.
`A1^3` and `A1^4` at `--maxdim 20` or higher are the most promising targets —
see below.

## Why `SL2^k` is the natural place to look

Section 3.4 of the paper conjectures that for `g` a product of copies of
`sl2`, the representations for which (2) holds are **exactly** the extended
quiver representations of trivial type. If the Hesselink-type class is
strictly larger there, the extra members are precisely the counterexamples,
and the gap between the two classes is the shape of the corrected conjecture.

## Status

Validated: reproduces Example 4.3 of the paper verbatim
(`X_V = 1 + q*chi_{V1 box V1} - q^2*chi_{V0 box V4}`), Example 4.2, and the
adjoint representations of `A1`, `A2`, `B2`, `G2`. The Theorem 1.2 detector is
calibrated in both directions (`V_1`, `V_1+V_1`, `V_2`, and the `A2` adjoint
must hold; `V_3`, `V_4` must fail) and reports exactly that.

**The `SL2^k` question is settled in range, and the answer is informative.**
For `SL2^2` (`dim V <= 12`) and `SL2^3` (`dim V <= 10`), with at most 4
summands:

    Hesselink-type  =  trivial-type extended quivers      exactly
    18 / 18  and  84 / 84,  in both directions, no exceptions

So §3.4's conjecture holds in this range — and, since the trivial-type
extended quivers are precisely the class Krylov–Wang can already prove, a
counterexample to Conjecture 1.6 **cannot** live in `SL2^k` at these
dimensions. That is what motivated `genrep.py`.

Tested directly against Theorem 1.2, all holding:

| group | Hesselink-type representations tested |
|---|---|
| `SL2`, `SL2^2` | 20 (`dim V <= 8`) |
| `SL2^3` | partial (`dim V <= 10`) |
| `A2` | 6 (all of them, `dim V <= 12`) |
| `A1xA2` | 20 (all of them, `dim V <= 12`) |

Coregularity (necessary for cofreeness) holds for every Hesselink-type
representation found so far, so none refutes that half either.

**No counterexample yet.** Krylov states one exists, so the search range is
still too small. The untried directions, in order of promise:

1. **Higher rank and larger `dim V`** — `A3`, `A1^2xA2`, `A2xA2`, and
   `SL2^4`. `genrep.py` handles any product of type-A factors.
2. **Non-type-A groups** — `B2`, `B3`, `C3`, `G2`. The `X_V` and identity-(2)
   tests already work there (`census.py`); only the direct Theorem 1.2 test
   does not, because `genrep.py` builds fundamental representations as
   exterior powers of the standard representation, which is type-A specific.
   Adding `B`/`C`/`D` means the standard representation plus a bilinear form;
   `G2` needs its 7-dimensional representation written out.
3. **Full cofreeness** — coregular plus equidimensional. Still untested, and
   still the one gap Proposition 2.1 cannot bridge.

### A bug worth knowing about

An earlier version of the Theorem 1.2 test only compared multiplicities at the
`lambda` occurring in `C[N_V]`. That silently skipped every `lambda` where
`M_q` predicts a nonzero multiplicity and `C[N_V]` has none — which is exactly
the asymmetry a counterexample would show. Both `mq.py` and `general12.py` now
range over all dominant weights either side can see. The fix makes the test
strictly stronger (it finds three mismatches for `V_3` where the old one found
one), and every result quoted above was produced after it.

## Files

| file | contents |
|---|---|
| `rootdata.py` | root systems (incl. non-simply-laced and products), Weyl groups, Freudenthal weight multiplicities |
| `hess.py` | the exact `X_V` computation and the Hesselink-type test |
| `series.py` | Hilbert series of `C[V]^G` and `C[V^{rho^vee}]^h`; the identity-(2) test; the coregularity test |
| `nullcone.py` | explicit `SL2^k` representations, invariants, and `C[N_V]` by linear algebra |
| `mq.py` | the q-Kostant partition function `P_q`, `M_q`, and the direct Theorem 1.2 test |
| `genrep.py` | explicit irreducibles for any product of type-A factors, via exterior powers and Cartan components |
| `general12.py` | the direct Theorem 1.2 test for products of type-A factors |
| `extquiver.py` | extended quivers of trivial type (§3.4), generated from the recursive definition |
| `census.py` | search driver via identity (2) |
| `sweep12.py` | search driver via Theorem 1.2 (decisive) |

Weights are stored in the basis of fundamental weights as tuples of
`Fraction`s, so `<lambda, alpha_i^vee>` is the `i`-th coordinate and
`s_i(lambda) = lambda - a_i alpha_i`. Heights `<lambda, rho^vee>` are computed
by passing to the `alpha`-basis via the inverse Cartan matrix and summing
coordinates. The Cartan matrices are checked at construction against their
symmetrisers (`d_i A_ij = d_j A_ji`), and representation dimensions have been
spot-checked against known values (`F4`: 26, `B3` spin: 8, `G2`: 7 and 14).
