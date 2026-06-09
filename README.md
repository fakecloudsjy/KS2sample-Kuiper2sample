# KS2sample and Kuiper2sample

Wolfram Language implementations of the two-sample Kolmogorov–Smirnov and Kuiper tests based on the methodology of the R package KSgeneral.

These functions provide exact permutation inference for continuous, discrete, and mixed distributions, including data with ties. In addition to reproducing the underlying methodology of KSgeneral, the Wolfram Language implementation leverages exact arithmetic, allowing p-values to be represented as exact rational numbers or evaluated to arbitrary precision.

The original KSgeneral package is available at:

* [KSgeneral GitHub repository](https://github.com/d-dimitrova/KSgeneral)
* [KSgeneral CRAN repository](https://CRAN.R-project.org/package=KSgeneral)

This repository provides Wolfram Language implementations of the two-sample procedures `KS2sample` and `Kuiper2sample` and serves as a companion repository for their submission to the Wolfram Function Repository.

## Functions

* `KS2sample` — computes weighted two-sample Kolmogorov–Smirnov statistics and their exact permutation p-values.
* `Kuiper2sample` — computes the two-sample Kuiper statistic and its exact permutation p-value.

## Method

Given two samples `X = {X₁, ..., Xₘ}` and `Y = {Y₁, ..., Yₙ}`, the observations are pooled into

`Z = {X₁, ..., Xₘ, Y₁, ..., Yₙ}`.

Under the null hypothesis that both samples arise from the same distribution, the permutation distribution is obtained by considering all `Binomial[m+n, m]` reallocations of the pooled sample into two groups of sizes `m` and `n`. The p-value is the proportion of reallocations whose test statistic is at least as extreme as the observed value.

## Features

* Exact permutation inference
* Exact rational p-values and arbitrary-precision numerical output
* Supports continuous, discrete, and mixed distributions
* Supports ties
* Weighted Kolmogorov–Smirnov statistics
* User-defined weight functions
* Pure Wolfram Language implementation

## Repository Contents

* `KS2sample.nb` — documentation notebook for the weighted two-sample KS test.
* `Kuiper2sample.nb` — documentation notebook for the two-sample Kuiper test.

## Function Repository

These functions are intended for submission to the Wolfram Function Repository.

## References

Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026). Efficient Exact Calculation of p-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests. To appear in Journal of Statistical Computation and Simulation.
