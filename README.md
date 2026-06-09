# KS2sample and Kuiper2sample

Wolfram Language implementations of the exact two-sample Kolmogorov–Smirnov
and Kuiper tests, companion to the R package
[KSgeneral](https://CRAN.R-project.org/package=KSgeneral) and to the paper:

> Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026). Efficient Exact
> Calculation of p-values of the Two-sample Kolmogorov–Smirnov and Kuiper
> Tests. *To appear in Journal of Statistical Computation and Simulation.*

The functions provide exact permutation p-values for continuous, discrete, and
mixed distributions, including data with ties. Because the Wolfram Language
performs exact rational arithmetic natively, p-values are returned as exact
rational numbers and can be evaluated to arbitrary numerical precision — a
property not available in the R/C++ implementation.

Both functions are submitted to the
[Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).

---

## Functions

### `KS2sample`

Computes the exact p-value of the (weighted) two-sample Kolmogorov–Smirnov
test.

**Data interface** — accepts raw samples and returns `{statistic, p-value}`:

```mathematica
KS2sample[data1, data2]
KS2sample[data1, data2, alternative]
KS2sample[data1, data2, alternative, weight]
```

**Low-level interface** — accepts precomputed inputs and returns the p-value:

```mathematica
KS2sample[m, n, q]
KS2sample[m, n, alternative, q]
KS2sample[m, n, alternative, multiplicity, q]
KS2sample[m, n, alternative, multiplicity, q, wvec]
```

| Argument | Description |
|---|---|
| `data1`, `data2` | Numeric lists; the two samples. May contain ties. |
| `alternative` | `"TwoSided"` (default), `"Greater"`, or `"Less"`. |
| `weight` | `0` (unweighted, default); `ν ∈ (0,1]` for W(t) = 1/[t(1−t)]^ν; or a user-supplied pure function. |
| `m`, `n` | Sample sizes. |
| `multiplicity` | Integer list of length k: count of each distinct pooled-sample value. `Total[multiplicity] == m+n`. |
| `q` | Observed value of the KS statistic. |
| `wvec` | Weight vector of length m+n−1: `wvec[[i]] = W(i/(m+n))`. |

**Examples:**

```mathematica
(* Unweighted two-sided test *)
KS2sample[data1, data2]
(* -> {0.233333, 17/143} *)

(* One-sided test with Anderson–Darling weight (ν = 1/2) *)
KS2sample[data1, data2, "Greater", 1/2]

(* Simplest low-level call: no ties, two-sided *)
KS2sample[5, 6, 3/10]
(* -> 69/77 *)

(* Low-level call with tied observations *)
KS2sample[120, 150, "TwoSided", {80, 70, 40, 80}, 1/10]
```

---

### `Kuiper2sample`

Computes the exact p-value of the two-sample Kuiper test.

**Data interface:**

```mathematica
Kuiper2sample[data1, data2]
```

**Low-level interface:**

```mathematica
Kuiper2sample[m, n, q]
Kuiper2sample[m, n, multiplicity, q]
```

The Kuiper statistic V = Δ⁺ + Δ⁻ is the sum of the supremum and infimum of
the empirical CDF difference. Its invariance to cyclic shifts makes it
particularly suited to circular and seasonal data.

**Examples:**

```mathematica
(* Data interface *)
Kuiper2sample[data1, data2]

(* Simplest low-level call: no ties *)
Kuiper2sample[5, 6, 1/2]
(* -> 17/21 *)

(* Low-level call with tied observations *)
Kuiper2sample[120, 150, {80, 70, 40, 80}, 11/60]
```

---

## Method

Given two samples X = {X₁, …, Xₘ} and Y = {Y₁, …, Yₙ}, the observations
are pooled. Under the null hypothesis that both samples arise from the same
distribution, the permutation distribution considers all Binomial[m+n, m]
reallocations of the pooled sample into groups of sizes m and n. The p-value
is the proportion of reallocations for which the test statistic is at least as
extreme as the observed value.

The p-value is computed via the recurrence of Nikiforov (1994), extended by
Dimitrova, Jia, Kaishev (2026) to support arbitrary weight functions and
tied observations. All arithmetic is exact; no floating-point tolerance
parameter is required.

The recurrence runs in O(mn) time for `KS2sample`. For `Kuiper2sample` the
complexity is O(mn · LCM[m,n]): when m = n this is O(m³); when m and n are
coprime (e.g. m = n+1), LCM[m,n] = mn and complexity is O(m²n²).

By Theorem 1 (KS) and Theorem 3 (Kuiper) of Dimitrova, Jia, Kaishev (2026),
the permutation p-value is asymptotically consistent with the unconditional
p-value for arbitrary underlying distributions F and G.

---

## Features

| | `KS2sample` | `Kuiper2sample` |
|---|---|---|
| Exact rational p-values | ✓ | ✓ |
| Continuous data | ✓ | ✓ |
| Discrete data / ties | ✓ | ✓ |
| Unequal sample sizes | ✓ | ✓ |
| Weighted statistic | ✓ | — |
| User-defined weight function | ✓ | — |
| One-sided alternatives | ✓ | — |
| Circular data | — | ✓ |

---

## Repository Contents

| File | Description |
|---|---|
| `KS2sample.nb` | Wolfram Function Repository definition notebook for `KS2sample`. |
| `Kuiper2sample.nb` | Wolfram Function Repository definition notebook for `Kuiper2sample`. |
| `README.md` | This file. |

---

## Installation

Once published to the Wolfram Function Repository, the functions can be loaded
with:

```mathematica
ResourceFunction["KS2sample"][data1, data2]
ResourceFunction["Kuiper2sample"][data1, data2]
```

Until publication, load the definition notebooks directly in Mathematica.

---

## Relation to KSgeneral

This repository is a Wolfram Language companion to the R package KSgeneral.
The algorithms are identical; the implementations differ in two respects:

- **Exact arithmetic.** The Wolfram Language performs rational arithmetic
  natively. The R/C++ implementation uses floating-point arithmetic and
  requires a tolerance parameter `tol` to handle boundary cases; the
  Mathematica implementation instead uses `Ceiling[x] - 1` for the upper
  column bound, which enforces strict inequality exactly.

- **Output.** The data interfaces return `{statistic, p-value}` as a list.
  The low-level interfaces return the p-value only.

| R function | Mathematica equivalent |
|---|---|
| `KS2sample` | `KS2sample[data1, data2, …]` |
| `KS2sample_Rcpp` / `KS2sample_c_Rcpp` | `KS2sample[m, n, alternative, multiplicity, q, wvec]` |
| `Kuiper2sample` | `Kuiper2sample[data1, data2]` |
| `Kuiper2sample_Rcpp` / `Kuiper2sample_c_Rcpp` | `Kuiper2sample[m, n, multiplicity, q]` |

KSgeneral is available at:
- [KSgeneral on CRAN](https://CRAN.R-project.org/package=KSgeneral)
- [KSgeneral on GitHub](https://github.com/d-dimitrova/KSgeneral)

---

## References

Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026). Efficient Exact Calculation
of p-values of the Two-sample Kolmogorov–Smirnov and Kuiper Tests. *To appear
in Journal of Statistical Computation and Simulation.*

Nikiforov, A. M. (1994). Algorithm AS 288: Exact Smirnov Two-Sample Tests for
Arbitrary Distributions. *Journal of the Royal Statistical Society, Series C*,
43(1), 265–270.

Kuiper, N. H. (1960). Tests concerning random points on a circle.
*Proceedings Koninklijke Nederlandse Akademie van Wetenschappen A*, 63, 38–47.

Maag, U. R., Stephens, M. A. (1968). The V_NM two-sample test. *Annals of
Mathematical Statistics*, 39(3), 923–935.
