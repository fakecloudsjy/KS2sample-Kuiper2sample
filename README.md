# KS2sample and Kuiper2sample

Wolfram Language implementations of exact two-sample Kolmogorov–Smirnov and Kuiper permutation tests for continuous, discrete, and mixed distributions.

## Functions

### KS2sample

Computes exact permutation p-values for weighted two-sample Kolmogorov–Smirnov statistics, including:

* Two-sided KS statistic
* One-sided greater alternative
* One-sided less alternative
* User-defined weight functions
* Exact finite-sample permutation inference

### Kuiper2sample

Computes exact permutation p-values for the two-sample Kuiper test, a circularly invariant alternative to the Kolmogorov–Smirnov test that is often more sensitive to differences in both location and spread.

## Method

## Method

Given two samples `X = {X₁, ..., Xₘ}` and `Y = {Y₁, ..., Yₙ}`, the tests first form the pooled sample

`Z = {X₁, ..., Xₘ, Y₁, ..., Yₙ}`.

The reported p-value is an exact permutation p-value obtained by considering all possible reallocations of the pooled sample into two groups of sizes `m` and `n`. Equivalently, the permutation distribution is computed over all `Binomial[m+n, m]` possible partitions of the pooled sample.

## Features

* Exact permutation p-values
* Supports ties and discrete data
* Continuous, discrete, and mixed distributions
* Weighted KS statistics
* User-specified weight functions
* Pure Wolfram Language implementation

## Repository Contents

* `KS2sample.nb` — documentation notebook for the weighted two-sample KS test.
* `Kuiper2sample.nb` — documentation notebook for the two-sample Kuiper test.

## Function Repository

These functions are intended for submission to the Wolfram Function Repository.

## References

* Kolmogorov (1933)
* Smirnov (1948)
* Kuiper (1960)
* Anderson and Darling (1952)
* Dimitrova, Jia, and Kaishev (2026)
