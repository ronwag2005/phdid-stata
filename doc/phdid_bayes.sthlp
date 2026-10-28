{smcl}
{* *! version 0.1.0  26oct2026}{...}
{vieweralsosee "phdid" "help phdid"}{...}
{vieweralsosee "phdid l0" "help phdid_l0"}{...}
{vieweralsosee "phdid postestimation" "help phdid_postestimation"}{...}
{viewerjumpto "Syntax" "phdid_bayes##syntax"}{...}
{viewerjumpto "Description" "phdid_bayes##description"}{...}
{viewerjumpto "Options" "phdid_bayes##options"}{...}
{viewerjumpto "Remarks" "phdid_bayes##remarks"}{...}
{viewerjumpto "Examples" "phdid_bayes##examples"}{...}
{viewerjumpto "Stored results" "phdid_bayes##results"}{...}
{title:Title}

{phang}
{bf:phdid bayes} {hline 2} The Dirichlet process partial-homogeneity estimator


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:phdid bayes} [{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Prior}
{synopt:{opt alpha(#)}}concentration parameter; default {cmd:alpha(1)}{p_end}
{synopt:{opt mu0(#)}}base measure mean; default median of the effects{p_end}
{synopt:{opt s0sq(#)}}base measure variance; default (10 sd)^2{p_end}

{syntab:Sampler}
{synopt:{opt iter:s(#)}}iterations per chain; default {cmd:iters(4000)}{p_end}
{synopt:{opt burn:in(#)}}discarded iterations; default {cmd:burnin(1000)}{p_end}
{synopt:{opt thin(#)}}keep every {it:#}th draw; default {cmd:thin(1)}{p_end}
{synopt:{opt chains(#)}}number of chains; default {cmd:chains(1)}{p_end}
{synopt:{opt init(rule)}}{opt singletons}, {opt pooled} or
{opt dispersed}{p_end}
{synopt:{opt marg:inal(rule)}}{opt exact} (default) or {opt diagonal}{p_end}
{synopt:{opt seed(string)}}set the random-number seed first{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}credible level; default {cmd:level(95)}{p_end}
{synopt:{opt nococlust}}suppress the co-clustering table{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{cmd:phdid bayes} is {cmd:eclass}.  {helpb phdid_setup:phdid setup} must be run
first.


{marker description}{...}
{title:Description}

{pstd}
{cmd:phdid bayes} fits the Dirichlet process mixture over partitions of the
cohort-time cells by collapsed Gibbs sampling, Algorithm 3 of Neal (2000).
Each sweep reassigns one cell at a time, integrating the group effect out of
the assignment move, so the sampler moves over partitions rather than over
partitions and effects jointly.

{pstd}
Because it averages over partitions rather than selecting one, its intervals
include the uncertainty about the grouping itself.  In the paper's simulations
that is worth 0.79-0.81 coverage against 0.57-0.62 for the l0 estimator when
the partition is genuinely uncertain, at no cost when the effects are well
separated.  This is the recommended route for inference.

{pstd}
The marginal likelihood in the assignment moves uses the exact joint
covariance of the first-stage effects, through the Woodbury identity of
equation (30) of the paper.  That is what the {cmd:marginal(exact)} default
means, and it is what the paper reports.


{marker options}{...}
{title:Options}

{phang}
{opt alpha(#)} is the Dirichlet process concentration parameter.  Larger
values favour more groups: the prior expected number is
{it:alpha} ln(1 + {it:K}/{it:alpha}), which the command reports so that the
prior and posterior group counts can be compared directly.  Trace the
dependence with {cmd:phdid sensitivity, alpha()}.

{phang}
{opt mu0(#)} and {opt s0sq(#)} are the mean and variance of the normal base
measure.  The defaults are data-scaled -- the median of the first-stage
effects and (10 times their standard deviation) squared -- matching the
paper's applications.  This is mildly empirical-Bayes; supply both for a prior
that does not look at the data.

{phang}
{opt iters(#)}, {opt burnin(#)} and {opt thin(#)} control the run.  The
defaults are deliberately modest; for published results use at least
{cmd:iters(20000) burnin(2000)} and check {cmd:phdid rhat}.

{phang}
{opt chains(#)} runs several chains.  With {cmd:init(dispersed)} this is what
makes {cmd:phdid rhat} able to detect a chain stuck in one corner of partition
space, which is the failure mode that matters here.

{phang}
{opt init(rule)} sets the starting partition: {opt singletons} puts every cell
in its own group, {opt pooled} puts them all together, and {opt dispersed}
varies the start by chain.

{phang}
{opt marginal(rule)} chooses how the covariance enters the assignment moves.
{opt exact} uses the full joint covariance and is the default.  {opt diagonal}
uses the faster shortcut of the original replication scripts, which discards
the cross-cell terms; it is provided for comparison, and
{cmd:phdid covcheck} quantifies what it costs on your design.

{phang}
{opt level(#)} sets the credible level for the reported intervals, which are
posterior percentiles.

{phang}
{opt seed(string)} sets the random-number seed, so the run is reproducible.

{phang}
{opt nococlust} suppresses the co-clustering table.


{marker remarks}{...}
{title:Remarks}

{dlgtab:Co-clustering, not a partition}

{pstd}
The posterior does not concentrate on one partition and should not be reported
as if it did.  On the paper's first application the most probable grouping
carries under a quarter of the posterior mass.  The co-clustering matrix --
the posterior probability that two cells share a group -- is the summary that
survives: it says which pairs the data really do bind together.  See
{cmd:phdid graph coclust}, and {cmd:phdid partition} for a point summary when
one is needed.

{dlgtab:What e(V) is}

{pstd}
{cmd:e(V)} is the posterior covariance of the cell effects across draws.  It
already includes uncertainty about the partition, so {helpb lincom} and
{helpb nlcom} on {cmd:e(b)} give partition-averaged inference.  For the
standard estimands use {cmd:phdid aggregate}, which applies the aggregation to
each draw rather than to the posterior mean, and reports percentile
intervals.

{dlgtab:Checking the sampler}

{pstd}
Run several dispersed chains and check {cmd:phdid rhat}, which reports
split-chain Gelman-Rubin statistics and effective sample sizes for every cell
effect {it:and} for the number of groups.  The group count is the diagnostic
that matters: chains can agree on the effects while disagreeing about how many
groups there are.  On small designs, {cmd:phdid enumerate} computes the exact
posterior and the sampler can be checked against it directly.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. set seed 20260301}{p_end}
{phang2}{cmd:. phdid bayes, alpha(1) iters(20000) burnin(2000) chains(4) init(dispersed)}{p_end}
{phang2}{cmd:. phdid rhat}{p_end}
{phang2}{cmd:. phdid aggregate dynamic}{p_end}
{phang2}{cmd:. phdid graph coclust}{p_end}

{pstd}A prior that does not look at the data:{p_end}

{phang2}{cmd:. phdid bayes, alpha(1) mu0(0) s0sq(1)}{p_end}

{pstd}What the diagonal shortcut would have given:{p_end}

{phang2}{cmd:. phdid bayes, marginal(diagonal)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:phdid bayes} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(m_mean)}}posterior mean number of groups{p_end}
{synopt:{cmd:e(prior_m)}}prior mean number of groups{p_end}
{synopt:{cmd:e(alpha)}}concentration parameter{p_end}
{synopt:{cmd:e(mu0)}, {cmd:e(s0sq)}}base measure actually used{p_end}
{synopt:{cmd:e(ndraws)}}draws kept{p_end}
{synopt:{cmd:e(chains)}, {cmd:e(iters)}, {cmd:e(burnin)}}run settings{p_end}
{synopt:{cmd:e(K)}}number of cells{p_end}
{synopt:{cmd:e(level)}}credible level{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:phdid_bayes}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(marginal)}}{cmd:exact} or {cmd:diagonal}{p_end}
{synopt:{cmd:e(init)}}starting rule{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}posterior mean cell effects{p_end}
{synopt:{cmd:e(V)}}their posterior covariance{p_end}
{synopt:{cmd:e(cells)}}flexible estimate, posterior mean, sd and
interval{p_end}
{synopt:{cmd:e(coclust)}}posterior co-clustering probabilities{p_end}
{p2colreset}{...}

{pstd}
The draws themselves are left in Mata as {cmd:PHDID_TAUDRAWS},
{cmd:PHDID_ZDRAWS}, {cmd:PHDID_MDRAWS} and {cmd:PHDID_CHAIN}, which is what
{cmd:phdid aggregate}, {cmd:phdid rhat} and {cmd:phdid partition} read.
