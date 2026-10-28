# phdid for Stata 0.1.0

First release, implementing Arora and Wagle (2026), "A Bayesian Approach to
Partial Homogeneity in Staggered Difference-in-Difference"
(doi:10.2139/ssrn.7207083).

## Input

* `phdid setup` assembles the inputs from any first-stage command that posts
  `e(b)` and `e(V)`, or from matrices supplied directly. Lemma 1 of the paper
  guarantees that two first stages giving the same `(b, V)` give the same
  partition posterior; the certification script asserts it.
* Coefficient names are parsed for the conventions used by `csdid` and the
  other common commands, with `cells()` as the fallback.
* `phdid describe` reports what is in memory and warns on a diagonal
  covariance.

## Estimation

* `phdid bayes` fits the Dirichlet process mixture by collapsed Gibbs
  sampling (Neal 2000, Algorithm 3), with the exact collapsed marginal
  likelihood in the cluster-assignment moves by default.
* `phdid l0` computes the two-step l0 estimator of Appendix B, supporting BIC
  selection along the agglomeration path, the fixed-`lambda` stopping rule,
  and a fixed group count.
* `phdid fit` gives the restricted GLS fit under a partition supplied by the
  user. Because the partition is not selected on the data, its deviance is a
  valid chi-squared test of the restriction.

## Inference

* `phdid aggregate` applies overall, event-study, by-cohort and by-calendar
  aggregators to each posterior draw, so reported intervals propagate
  uncertainty about the partition itself. After `phdid l0` or `phdid fit` the
  same aggregators are applied by the delta method.
* `phdid partition` summarises the posterior by a point partition under
  Wade-Ghahramani variation of information or Binder loss, and reports how
  little posterior mass that single partition carries.

## Specification testing

* `phdid test` reports a model-implied common-effect chi-squared test on the
  pooled GLS deviance, an excess-dispersion heterogeneity share computed
  against the exact covariance, and — when placebo estimates are supplied —
  the pre/post noise gauge and within-cell randomization test of section
  5.2.3. Each block is tagged with its provenance.

## Diagnostics

* `phdid enumerate` computes the exact partition posterior for small designs,
  the Appendix E benchmark, and is the oracle the sampler is certified
  against.
* `phdid rhat` reports split-chain Gelman-Rubin statistics and effective
  sample sizes across dispersed chains, for the cell effects and for the
  number of groups.
* `phdid covcheck` compares exact against diagonal handling of the covariance
  in the cluster-assignment moves.

## Regularisation paths

* `phdid sensitivity` traces the DP estimator across the concentration
  parameter and the l0 estimator across its penalty. `type(cells)` reports the
  path for every cohort-time effect separately rather than only for the
  aggregate; this is usually the more informative view, since the overall ATT
  is robust to over-pooling and a flat aggregate path can hide cells that move
  a great deal.
* `by(m)` walks the agglomeration path one merge at a time. This is not the
  same as varying `lambda`: the Appendix B stopping rule compares each merge
  cost against `lambda * |A| * |B|`, and that factor grows as groups absorb
  one another, so the effective per-merge threshold is not monotone in the
  merge order and some group counts are unreachable by any single `lambda`. On
  the paper's first application the thresholds run 0.079, 0.686, 0.397, 1.350,
  2.765, 3.571, so nothing selects five groups.

## Graphs

* `phdid graph coclust`, `groups` and `bic`, plus a graph from
  `phdid sensitivity`. All are built from `twoway` primitives with no
  dependency on user-written graphing packages. The co-clustering grid is
  drawn with `rbar` at one data unit per cell rather than with `contour`,
  which would interpolate between grid points and imply a continuum between
  cohort-time cells that does not exist.

## Notes on choices that move results

* The BIC penalty on the two-stage route needs a sample size that the theory
  does not fix. The default is `K`, which reproduces the paper's reported
  group counts; `nbic()` exposes the choice because it changes the selected
  partition.
* The default base measure is data-scaled (`mu0` = median of the effects,
  `s0sq` = (10 sd)^2), matching the paper's applications. This is mildly
  empirical-Bayes; supply both for a prior independent of the data.
* `phdid bayes, marginal(diagonal)` reproduces the faster shortcut used in the
  original replication scripts. The default is `exact`, matching what the
  paper reports.
