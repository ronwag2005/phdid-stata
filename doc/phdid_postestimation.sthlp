{smcl}
{* *! version 0.1.0  26oct2026}{...}
{vieweralsosee "phdid" "help phdid"}{...}
{vieweralsosee "phdid l0" "help phdid_l0"}{...}
{vieweralsosee "phdid bayes" "help phdid_bayes"}{...}
{viewerjumpto "Syntax" "phdid_postestimation##syntax"}{...}
{viewerjumpto "aggregate" "phdid_postestimation##aggregate"}{...}
{viewerjumpto "fit" "phdid_postestimation##fit"}{...}
{viewerjumpto "partition" "phdid_postestimation##partition"}{...}
{viewerjumpto "rhat" "phdid_postestimation##rhat"}{...}
{viewerjumpto "sensitivity" "phdid_postestimation##sensitivity"}{...}
{viewerjumpto "covcheck" "phdid_postestimation##covcheck"}{...}
{viewerjumpto "enumerate" "phdid_postestimation##enumerate"}{...}
{viewerjumpto "graph" "phdid_postestimation##graph"}{...}
{title:Title}

{phang}
{bf:phdid postestimation} {hline 2} Aggregation, diagnostics and graphs after
{cmd:phdid}


{marker syntax}{...}
{title:Syntax}

{synoptset 34 tabbed}{...}
{synopthdr:subcommand}
{synoptline}
{synopt:{helpb phdid_postestimation##aggregate:aggregate}}overall,
event-study, by-cohort and by-period estimands{p_end}
{synopt:{helpb phdid_postestimation##fit:fit}}restricted GLS under a partition
you supply{p_end}
{synopt:{helpb phdid_postestimation##partition:partition}}a point partition
from the posterior{p_end}
{synopt:{helpb phdid_postestimation##rhat:rhat}}sampler convergence{p_end}
{synopt:{helpb phdid_postestimation##sensitivity:sensitivity}}regularisation
paths{p_end}
{synopt:{helpb phdid_postestimation##covcheck:covcheck}}exact against diagonal
covariance{p_end}
{synopt:{helpb phdid_postestimation##enumerate:enumerate}}the exact partition
posterior{p_end}
{synopt:{helpb phdid_postestimation##graph:graph}}co-clustering, grouped
effects, selection path{p_end}
{synoptline}
{p2colreset}{...}


{marker aggregate}{...}
{title:phdid aggregate}

{p 8 15 2}
{cmd:phdid aggregate} [{it:estimand}] [{cmd:,} {opt l:evel(#)} {opt post}]

{pstd}
{it:estimand} is {opt overall} (the default), {opt dynamic} (by event time),
{opt group} (by cohort) or {opt calendar} (by period).

{pstd}
Every estimand the package reports is a fixed linear map of the cell effects.
After {cmd:phdid bayes} that map is applied to {it:each posterior draw}, so the
reported interval includes uncertainty about the grouping itself; the
intervals are posterior percentiles.  After {cmd:phdid l0} or
{cmd:phdid fit} it is applied to {cmd:e(b)} and {cmd:e(V)} by the delta
method, and the intervals condition on the partition.

{pstd}
Weights within each term come from {cmd:phdid setup} and are renormalised to
sum to one within the term.  {opt dynamic}, {opt group} and {opt calendar}
require that the cohort and period of each cell were recovered; see
{cmd:cells()} in {helpb phdid_setup:phdid setup}.

{phang}
{opt post} replaces {cmd:e(b)} and {cmd:e(V)} with the aggregates, so that
{helpb lincom} and {helpb test} operate on the terms shown.  This discards the
cell-level fit; re-run the estimation command to get it back.

{pstd}
Returns {cmd:r(table)} (estimate, standard error, lower, upper),
{cmd:r(A)} (the aggregation matrix), {cmd:r(b)}, {cmd:r(V)} and
{cmd:r(nterms)}.


{marker fit}{...}
{title:phdid fit}

{p 8 15 2}
{cmd:phdid fit}{cmd:,} {cmdab:part:ition(}{it:numlist}{cmd:)}
[{opt l:evel(#)} {opt notable}]

{pstd}
Fits restricted GLS under a partition supplied as one group number per cell,
in the order the cells appear in {cmd:phdid describe}.  The group numbers are
names, so {cmd:partition(1 1 2)} and {cmd:partition(3 3 1)} give identical
fits.

{pstd}
This is the route for a grouping that comes from theory, from an earlier
sample, or from {cmd:phdid partition}.  Because the partition was not selected
on these data, the reported deviance {it:is} a valid chi-squared test of the
restriction against the flexible fit -- unlike the same statistic after
{cmd:phdid l0}.

{pstd}
{cmd:eclass}.  Posts {cmd:e(b)}, {cmd:e(V)}, {cmd:e(cells)}, {cmd:e(phi)},
{cmd:e(partition)}, and the scalars {cmd:e(m)}, {cmd:e(deviance)},
{cmd:e(df)}, {cmd:e(p)}, {cmd:e(att)}, {cmd:e(att_se)}, {cmd:e(K)}.


{marker partition}{...}
{title:phdid partition}

{p 8 15 2}
{cmd:phdid partition} [{cmd:,} {cmd:loss(}{opt vi}|{opt binder}{cmd:)}
{opt post}]

{pstd}
Summarises the posterior over partitions by a single grouping, minimising an
expected loss over the partitions the sampler visited.  {opt vi}, the default,
uses the variation of information of Wade and Ghahramani (2018), which
respects the lattice structure of partition space.  {opt binder} uses Binder's
loss, the L1 distance to the co-clustering matrix.

{pstd}
The posterior mode is deliberately not offered: with {it:K} cells no single
partition holds much mass, and the command reports how little the chosen one
holds so that this is visible.  A point partition is a summary, not a finding;
the interval from {cmd:phdid aggregate} remains the one to report.

{phang}
{opt post} posts the restricted GLS fit under the chosen partition.  Those are
conditional standard errors.

{pstd}
Returns {cmd:r(partition)}, {cmd:r(coclust)}, {cmd:r(b)}, {cmd:r(V)},
{cmd:r(m)}, {cmd:r(loss_value)}, {cmd:r(post_prob)} and
{cmd:r(partition_str)}.  Requires {cmd:phdid bayes}.


{marker rhat}{...}
{title:phdid rhat}

{p 8 15 2}
{cmd:phdid rhat} [{cmd:,} {opt maxl:ag(#)}]

{pstd}
Split-chain Gelman-Rubin statistics and effective sample sizes for every cell
effect and for the number of groups, the latter computed by Geyer's initial
positive sequence truncated at {opt maxlag()} (default 200).

{pstd}
The number of groups is the diagnostic that matters here.  The partition label
is not a continuous parameter, so a chain stuck in one corner of partition
space can still look settled on the cell effects.  Run at least four chains
with {cmd:init(dispersed)} for this to have power.

{pstd}
Returns {cmd:r(table)}, {cmd:r(rhat_max)}, {cmd:r(ess_min)} and
{cmd:r(nchains)}.  Requires {cmd:phdid bayes}.


{marker sensitivity}{...}
{title:phdid sensitivity}

{p 8 15 2}
{cmd:phdid sensitivity} [{cmd:,} {opt alpha(numlist)} {opt lambda(numlist)}
{cmd:by(}{opt lambda}|{opt m}{cmd:)} {cmd:type(}{opt overall}|{opt cells}{cmd:)}
{opt l:evel(#)} {opt saving(filename)} {opt nograph}
{opt graphopts(string)} {it:sampler_options}]

{pstd}
Re-estimates across the tuning parameter and reports the path.  An explicit
{opt alpha()} traces the Dirichlet process estimator; an explicit
{opt lambda()}, or {cmd:by(m)}, traces the l0 estimator; with neither, the
command follows whichever estimator was run last.  {cmd:by(m)} walks the
agglomeration path one merge at a time, which is the way to reach group counts
that no single {it:lambda} selects.

{pstd}
{cmd:type(cells)} reports the path for every cohort-time effect separately
rather than for the aggregate, and is usually the more informative view: the
overall ATT is robust to over-pooling, so a flat aggregate path can hide cells
that move a great deal.

{pstd}
A graph is drawn unless {opt nograph} is specified; {opt graphopts()} passes
options to it, and {opt saving()} saves the path as a dataset.  The sampler
options {opt iters()}, {opt burnin()}, {opt thin()}, {opt chains()},
{opt marginal()}, {opt mu0()}, {opt s0sq()}, {opt init()} and {opt seed()} are
accepted and apply at every point on an {opt alpha()} path.

{pstd}
Returns {cmd:r(path)}, {cmd:r(npoints)}, {cmd:r(parameter)} and
{cmd:r(path_type)}; the l0 path also returns {cmd:r(thresholds)}, the
per-merge {it:lambda} implied by the agglomeration path.


{marker covcheck}{...}
{title:phdid covcheck}

{p 8 15 2}
{cmd:phdid covcheck} [{cmd:,} {opt alpha(#)} {opt mu0(#)} {opt s0sq(#)}
{opt iter:s(#)} {opt burn:in(#)} {opt thin(#)} {opt chains(#)} {opt init(rule)}
{opt seed(string)}]

{pstd}
Runs the sampler twice on the same data, once with the exact collapsed
marginal likelihood and once with the diagonal shortcut, and reports how far
apart they land: the ratio of posterior variances of the overall effect, and
the largest discrepancy in any co-clustering probability.

{pstd}
This is Remark 1 of the paper made operational.  Cohort-time cells share units
and share the two-way fixed effects, so the shortcut is not a harmless
speed-up.  In the paper's design it understates the posterior variance of
aggregates by about 3.5 times and can move a single co-clustering probability
by 0.39.

{pstd}
Returns {cmd:r(var_ratio)}, {cmd:r(cc_maxdiff)}, {cmd:r(cells)},
{cmd:r(cc_exact)}, {cmd:r(cc_diag)} and the per-run means.


{marker enumerate}{...}
{title:phdid enumerate}

{p 8 15 2}
{cmd:phdid enumerate} [{cmd:,} {opt alpha(#)} {opt mu0(#)} {opt s0sq(#)}
{opt maxc:ells(#)} {opt l:evel(#)} {opt nococlust}]

{pstd}
Computes the exact partition posterior by enumerating every partition of the
cells, the Appendix E benchmark.  The number of partitions is the Bell number,
which grows faster than exponentially -- 877 for 7 cells, 4,213,597 for 15 --
so {opt maxcells()} (default 11) refuses designs that would not finish.

{pstd}
This is the deterministic oracle the Gibbs sampler is validated against, and
the reason a different random-number generator in a different language does
not matter: both implementations compute this posterior exactly, and the
sampler is required only to reproduce it to within Monte Carlo error.

{pstd}
Returns {cmd:r(n_partitions)}, {cmd:r(m_mean)}, {cmd:r(att)},
{cmd:r(att_lo)}, {cmd:r(att_hi)}, {cmd:r(coclust)} and {cmd:r(tau)}.


{marker graph}{...}
{title:phdid graph}

{p 8 15 2}
{cmd:phdid graph} [{it:kind}] [{cmd:,} {it:twoway_options}]

{pstd}
{it:kind} is one of:

{phang2}
{opt coclust} draws the posterior co-clustering matrix as a shaded grid, one
cell per pair, with the probabilities printed.  Each cell is drawn as a bar
one data unit wide, so the grid is exact; contour plots were not used because
they interpolate between grid points, implying a continuum between cohort-time
cells that does not exist.  {opt novalues} suppresses the printed numbers.

{phang2}
{opt groups} plots the flexible first-stage estimates against the grouped or
posterior-mean estimates, cell by cell, with intervals.  The point of the
estimator is visible here: within a recovered group the grouped estimates
collapse onto one value with a shorter interval.

{phang2}
{opt bic} draws the l0 selection path, deviance and BIC against the number of
groups, with the selected {it:m} marked.

{pstd}
With no {it:kind}, {opt coclust} is drawn after {cmd:phdid bayes} and
{opt groups} otherwise.  All three accept the usual {help twoway_options}.
