# KS2Sample and Kuiper2Sample

Wolfram Language implementations of the exact two-sample Kolmogorov–Smirnov
and Kuiper tests, companion to the R package
[KSgeneral](https://CRAN.R-project.org/package=KSgeneral) and to the paper:

> Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026). Efficient Exact
> Calculation of p-values of the Two-sample Kolmogorov–Smirnov and Kuiper
> Tests. *Journal of Statistical Computation and Simulation*, 1–19.
> doi:10.1080/00949655.2026.2721410

The functions provide exact permutation p-values for continuous, discrete, and
mixed distributions, including data with ties. Because the Wolfram Language
performs exact rational arithmetic natively, p-values are returned as exact
rational numbers and can be evaluated to arbitrary numerical precision — a
property not available in the R/C++ implementation in `KSgeneral`, which uses
floating-point arithmetic.

Both functions are published on the Wolfram Function Repository:

- [`KS2Sample`](https://resources.wolframcloud.com/FunctionRepository/resources/KS2Sample/)
- [`Kuiper2Sample`](https://resources.wolframcloud.com/FunctionRepository/resources/Kuiper2Sample/)

---

## Repository Contents

| Path | Description |
|---|---|
| `KS2Sample/KS2Sample.nb` | Wolfram Function Repository definition notebook for `KS2Sample`. |
| `KS2Sample/KS2Sample.wl` | Plain-text Wolfram Language implementation of `KS2Sample`. |
| `Kuiper2Sample/Kuiper2Sample.nb` | Wolfram Function Repository definition notebook for `Kuiper2Sample`. |
| `Kuiper2Sample/Kuiper2Sample.wl` | Plain-text Wolfram Language implementation of `Kuiper2Sample`. |
| `CITATION.cff` | Citation metadata for GitHub and citation-management tools. |
| `LICENSE` | MIT License for the software in this repository. |
| `README.md` | Project documentation. |

---

## Functions

### `KS2Sample`

Computes exact permutation p-values for the (weighted) two-sample
Kolmogorov–Smirnov test. Supports continuous, discrete, and mixed
distributions, including data with ties.

**Data interface** — accepts raw samples, returns `{statistic, p-value}`:

```mathematica
KS2Sample[data1, data2]
KS2Sample[data1, data2, "Alternative" -> alt]
KS2Sample[data1, data2, "Weight" -> w]
KS2Sample[data1, data2, "Alternative" -> alt, "Weight" -> w]
```

**Low-level interface** — accepts precomputed inputs, returns the p-value:

```mathematica
KS2Sample[m, n, q]
KS2Sample[m, n, multiplicity, q]
KS2Sample[m, n, multiplicity, q, "Alternative" -> alt]
KS2Sample[m, n, multiplicity, q, "Weight" -> w]
```

| Argument / Option | Description |
|---|---|
| `data1`, `data2` | Numeric lists; the two samples. May contain ties. |
| `"Alternative"` | `"TwoSided"` (default), `"Greater"`, or `"Less"`. |
| `"Weight"` | `0` (unweighted, default); `ν ∈ (0,1]` for W(t) = 1/[t(1−t)]^ν; or a user-supplied pure function. |
| `m`, `n` | Sample sizes. |
| `multiplicity` | Integer list: count of each distinct pooled-sample value. `Total[multiplicity] == m+n`. |
| `q` | Observed value of the KS statistic. Exact or approximate. |

**Examples:**

```mathematica
(* Unweighted two-sided test *)
KS2Sample[data1, data2]

(* One-sided test with Anderson–Darling weight (ν = 1/2) *)
KS2Sample[data1, data2, "Alternative" -> "Greater", "Weight" -> 1/2]

(* Low-level call: no ties, two-sided *)
KS2Sample[5, 6, 3/10]
(* -> 69/77 *)

(* Low-level call with tied observations *)
KS2Sample[120, 150, {80, 70, 40, 80}, 1/10]
```

---

### `Kuiper2Sample`

Computes exact permutation p-values for the two-sample Kuiper test. The
Kuiper statistic V = Δ⁺ + Δ⁻ is the sum of the two one-sided suprema of the
empirical CDF difference. Its invariance under cyclic shifts makes it
particularly suited to circular and seasonal data. Supports continuous,
discrete, and mixed distributions, including data with ties.

**Data interface:**

```mathematica
Kuiper2Sample[data1, data2]
```

**Low-level interface:**

```mathematica
Kuiper2Sample[m, n, q]
Kuiper2Sample[m, n, multiplicity, q]
```

**Examples:**

```mathematica
(* Data interface *)
Kuiper2Sample[data1, data2]

(* Low-level call: no ties *)
Kuiper2Sample[5, 6, 1/2]
(* -> 17/21 *)

(* Low-level call with tied observations *)
Kuiper2Sample[120, 150, {80, 70, 40, 80}, 11/60]
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
Dimitrova, Jia, Kaishev (2026) to support arbitrary weight functions and tied
observations. All arithmetic is exact; no floating-point tolerance parameter
is required.

The recurrence runs in O(mn) time for `KS2Sample`. For `Kuiper2Sample` the
complexity is O(mn · LCM[m,n]): when m = n this is O(m³); when m and n are
coprime (e.g. m = n+1), LCM[m,n] = mn and complexity is O(m²n²).

By Theorem 2.1 (KS) and Theorem 2.6 (Kuiper) of Dimitrova, Jia, Kaishev (2026),
the permutation p-value is asymptotically consistent with the unconditional
p-value for arbitrary underlying distributions F and G.

---

## Features

| | `KS2Sample` | `Kuiper2Sample` |
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

## Installation

Both functions are available directly from the Wolfram Function Repository:

- [`KS2Sample`](https://resources.wolframcloud.com/FunctionRepository/resources/KS2Sample/)
- [`Kuiper2Sample`](https://resources.wolframcloud.com/FunctionRepository/resources/Kuiper2Sample/)

In a Mathematica or Wolfram Language session, call them via `ResourceFunction`:

```mathematica
ResourceFunction["KS2Sample"][data1, data2]
ResourceFunction["Kuiper2Sample"][data1, data2]
```

No installation is required; `ResourceFunction` fetches and caches the function
automatically on first use.

### Wolfram Language Source

The `.wl` files in this repository provide the implementations as plain-text
Wolfram Language source, suitable for inspection, version control, or direct
use without a network connection:

```mathematica
Get["KS2Sample/KS2Sample.wl"]
Get["Kuiper2Sample/Kuiper2Sample.wl"]
```

---

## Relation to KSgeneral

This repository is a Wolfram Language companion to the R package KSgeneral.
The algorithms are identical; the implementations differ in two respects:

- **Exact arithmetic.** The Wolfram Language performs rational arithmetic
  natively. The R/C++ implementation uses floating-point arithmetic and
  requires a tolerance parameter `tol`; the Mathematica implementation uses
  `Ceiling[x] - 1` for the upper column bound, enforcing strict inequality
  exactly without any tolerance.

- **Interface.** The data interfaces return `{statistic, p-value}` as a list.
  The low-level interfaces return the p-value only. `"Alternative"` and
  `"Weight"` are named options rather than positional arguments.

| R function | Mathematica equivalent |
|---|---|
| `KS2sample` | `KS2Sample[data1, data2, …]` |
| `KS2sample_Rcpp` / `KS2sample_c_Rcpp` | `KS2Sample[m, n, multiplicity, q, …]` |
| `Kuiper2sample` | `Kuiper2Sample[data1, data2]` |
| `Kuiper2sample_Rcpp` / `Kuiper2sample_c_Rcpp` | `Kuiper2Sample[m, n, multiplicity, q]` |

KSgeneral is available at:
- [KSgeneral on CRAN](https://CRAN.R-project.org/package=KSgeneral)
- [KSgeneral on GitHub](https://github.com/d-dimitrova/KSgeneral)

---

## References

Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026). Efficient Exact Calculation
of p-values of the Two-sample Kolmogorov–Smirnov and Kuiper Tests. *Journal
of Statistical Computation and Simulation*, 1–19.
doi:10.1080/00949655.2026.2721410

Nikiforov, A. M. (1994). Algorithm AS 288: Exact Smirnov Two-Sample Tests for
Arbitrary Distributions. *Journal of the Royal Statistical Society, Series C*,
43(1), 265–270.

Finner, H., Gontscharuk, V. (2018). Two-sample Kolmogorov–Smirnov-type tests
revisited: Old and new tests in terms of local levels. *The Annals of
Statistics*, 46(6A), 3014–3037. doi:10.1214/17-AOS1647

Büning, H. (2001). Kolmogorov–Smirnov- and Cramér-von Mises Type Two-sample
Tests With Various Weight Functions. *Communications in Statistics —
Simulation and Computation*, 30(4), 847–865. doi:10.1081/SAC-100107780

Canner, P. L. (1975). A Simulation Study of One- and Two-Sample
Kolmogorov–Smirnov Statistics with a Particular Weight Function. *Journal of
the American Statistical Association*, 70(349), 209–211.

Kuiper, N. H. (1960). Tests concerning random points on a circle.
*Proceedings Koninklijke Nederlandse Akademie van Wetenschappen A*, 63, 38–47.

Maag, U. R., Stephens, M. A. (1968). The V_NM two-sample test. *Annals of
Mathematical Statistics*, 39(3), 923–935.

Hirakawa, K. (1973). The two-sample Kuiper test. *TRU Mathematics*, 9, 99–118.
