# phdid for Stata

**Partial homogeneity in staggered difference-in-differences.**

In a staggered DiD design the treatment effect is not a single number but a
vector of cohort-time effects (CATTs), one per cohort-time cell. Estimating
every one separately is unbiased but inefficient when some are in fact equal;
pooling them all into a single TWFE coefficient is efficient but biased
whenever the heterogeneity is genuine.

`phdid` treats the choice between these extremes as a **partition-selection
problem** on the cells, and implements the two estimators of Arora and Wagle
(2026) that solve it, along with the specification tests and diagnostics
needed to know whether either should be used at all.

This is the Stata implementation. An R package of the same name, with the same
estimators, is on CRAN. The two have been certified against each other: see
[Validation](#validation).

---

## Installation

```stata
net install phdid, from("https://raw.githubusercontent.com/ronwag2005/phdid-stata/main")
```

Requires Stata 17 or later. There are no other dependencies: the numerical
core is Mata, compiled into `lphdid.mlib`, with no external plugins.

To run the examples you will also want a heterogeneity-robust first-stage
command such as [`csdid`](https://github.com/friosavila/stpackages).

## The idea in sixty seconds

With `G` treated cohorts and `T` periods there are `K` cohort-time cells. The
true effect vector exhibits *partial homogeneity* if there is a partition of
those cells into `m < K` groups within which the effects are equal. Under the
true partition the grouped estimator is both unbiased and strictly more
efficient than the flexible one — in a balanced design, exactly `|C|` times
more efficient for every cell in a group of size `|C|`.

The partition is unknown, so the package recovers it two ways:

| | `phdid l0` | `phdid bayes` |
|---|---|---|
| **Returns** | one partition | a posterior over partitions |
| **Selects by** | BIC on the agglomeration path, or a fixed `lambda` | Dirichlet process prior, collapsed Gibbs |
| **Intervals** | condition on the selected partition | marginalize over the partition |
| **Coverage when the partition is uncertain** | 0.57–0.62 | 0.79–0.81 |
| **Coverage when effects are separated** | 0.95 | 0.94 |

The two are not rivals: `phdid l0` is the fixed-variance MAP of the same
Bayesian model, differing only in the partition prior. Use `phdid l0` when you
want one interpretable grouping to report and `phdid bayes` for inference.

## Quick start

```stata
webuse set https://friosavila.github.io/playingwithstata/drdid/
webuse mpdta, clear

* Any heterogeneity-robust first stage will do.
csdid lemp, ivar(countyreal) time(year) gvar(first_treat) method(dripw) notyet

phdid setup          // keeps post-treatment cells + the exact joint covariance

* 1. Is there heterogeneity to recover at all?
phdid test
//  chi2(6) = 25.29,  Prob > chi2 = 0.0003
//  heterogeneity share: 80%
//  Reading: recoverable heterogeneity.

* 2. One partition to report.
phdid l0
//  Groups (m) : 4  (selected by BIC, n = 7 in the penalty)
phdid aggregate
//  overall  -0.03883  [-0.06143, -0.01623]

* 3. Inference that accounts for not knowing the partition.
set seed 20260301
phdid bayes, alpha(1) iters(20000) burnin(2000) chains(4) init(dispersed)
phdid rhat
phdid aggregate
//  overall  -0.02382  [-0.04763, -0.00002]

phdid graph coclust  // which groupings are firm
```

## Why the Stata implementation is a faithful one

Lemma 1 of the paper says everything the estimators need enters only through
`R'Σ⁻¹τ` and `R'Σ⁻¹R`. That maps exactly onto Stata's universal `e(b)`/`e(V)`
convention, so `phdid setup` works after *any* first-stage command that posts
coefficients and a covariance — `csdid`, `jwdid`, `did_imputation`, or a
matrix you built yourself. Nothing is reimplemented from micro data, and
nothing is assumed about how the first stage was computed.

## Two things the package will keep telling you

**The covariance is not diagonal.** Cohort-time cells share units and share
the unit and time fixed effects, so their estimates are correlated. Using the
exact joint covariance is what makes the reported intervals honest: treating
the cells as independent understates the posterior variance of aggregates by
about 3.5× in the paper's design and can misstate individual co-clustering
probabilities by up to 0.39. `phdid setup` warns when handed a diagonal
covariance, and `phdid covcheck` quantifies the difference on your design.

**A partition is only worth reporting when the effects are separated enough to
be recovered.** Below that threshold both estimators return groupings that are
quantiles of noise. In the paper's second application the l0 estimator splits
138 event effects into five neat "bands" while a randomization test cannot
reject a single common effect. Run `phdid test` first; it says so plainly
rather than letting the output speak for itself.

## What's included

**Input** — `phdid setup` (from `e(b)`/`e(V)`, or from matrices you supply),
`phdid describe`.

**Specification testing** — `phdid test`, bundling a model-implied
common-effect chi-squared test, an excess-dispersion heterogeneity share
computed against the exact covariance, and the placebo gauge and within-cell
randomization test of section 5.2.3.

**Estimation** — `phdid l0`, `phdid bayes`, `phdid fit` (restricted GLS under
a partition you supply, whose deviance *is* a valid test since the partition
was not selected).

**Inference** — `phdid aggregate` for overall, event-study, by-cohort and
by-calendar estimands, computed per posterior draw after `phdid bayes`;
`phdid partition` for a point partition under variation-of-information or
Binder loss.

**Diagnostics** — `phdid enumerate` for the exact posterior on small designs,
`phdid rhat` for multi-chain convergence, `phdid sensitivity` for the
regularisation paths (per-cell as well as aggregate), `phdid covcheck` for
exact vs. diagonal.

**Graphs** — `phdid graph coclust`, `groups`, `bic`, and a graph from
`phdid sensitivity`. All built from `twoway` primitives, with no dependency on
user-written graphing packages.

## Validation

Everything the estimators do is deterministic given `(τ, Σ)` except the Gibbs
sampler, and Stata's random-number generator is not R's. So the package is
certified in two layers:

- **The deterministic parts** — the restricted GLS fit, the agglomeration
  path, BIC selection, the specification test, and the exact partition
  posterior by enumeration — are certified against the CRAN R package,
  using 143 reference values written at 17 significant figures.
- **The Gibbs sampler** is certified against the exact partition posterior
  from full enumeration, which both languages compute deterministically and
  which is the sampler's own target. That is the right benchmark: reproducing
  R's draws would be meaningless, reproducing R's posterior is not.

```stata
do cert/phdid_cert.do
//  phdid certification: 165 checks, all passed.
```

The script also checks invariants that do not depend on R: relabelling a
partition does not change the fit, aggregation weights are a partition of
unity, grouping never raises a cell's sampling variance, and a coarser
partition never fits better.

## Building from source

```stata
do build.do      // compiles src/*.mata into src/lphdid.mlib
```

The `.mata` sources ship for inspection and for certification; an installed
copy runs from the compiled library.

## Authorship and citation

The methods implemented here are joint work by **Parush Arora** and **Rohan
Wagle** (Department of Economics, Ashoka University). Both are copyright
holders under the MIT license. **Rohan Wagle** is the sole maintainer.

If you use `phdid`, please cite **both** the Ashoka University discussion
paper and the SSRN working paper:

> Arora, Parush and Rohan Wagle (2026). "A Bayesian Approach to Partial
> Homogeneity in Staggered Difference-in-Difference." Ashoka University
> Economics Discussion Paper 166.
> [Link](https://www.ashoka.edu.in/research/a-bayesian-approach-to-partial-homogeneity-in-staggered-difference-in-difference/)

> Arora, Parush and Rohan Wagle (2026). "A Bayesian Approach to Partial
> Homogeneity in Staggered Difference-in-Difference" (March 01, 2026). SSRN
> Working Paper 7207083. [SSRN](https://ssrn.com/abstract=7207083) |
> doi:[10.2139/ssrn.7207083](https://doi.org/10.2139/ssrn.7207083)

## References

Arora, P. and Wagle, R. (2026). A Bayesian Approach to Partial Homogeneity in
Staggered Difference-in-Difference. Ashoka University Economics Discussion
Paper 166.

Arora, P. and Wagle, R. (2026). A Bayesian Approach to Partial Homogeneity in
Staggered Difference-in-Difference (March 01, 2026). SSRN Working Paper
7207083. <https://doi.org/10.2139/ssrn.7207083>

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
Multiple Time Periods. *Journal of Econometrics* 225(2), 200–230.

Neal, R. M. (2000). Markov Chain Sampling Methods for Dirichlet Process
Mixture Models. *JCGS* 9(2), 249–265.

Wade, S. and Ghahramani, Z. (2018). Bayesian Cluster Analysis: Point
Estimation and Credible Balls. *Bayesian Analysis* 13(2), 559–626.

Wooldridge, J. M. (2025). Two-Way Fixed Effects, the Two-Way Mundlak
Regression, and Difference-in-Differences Estimators. *Empirical Economics*
69, 2545–2587.

## License

MIT © 2026 Rohan Wagle and Parush Arora. See [LICENSE](LICENSE).
