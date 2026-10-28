{smcl}
{* *! version 0.1.0  29oct2026}{...}
{vieweralsosee "[R] estimation" "help estcom"}{...}
{vieweralsosee "phdid setup" "help phdid_setup"}{...}
{vieweralsosee "phdid l0" "help phdid_l0"}{...}
{vieweralsosee "phdid bayes" "help phdid_bayes"}{...}
{vieweralsosee "phdid test" "help phdid_test"}{...}
{vieweralsosee "phdid postestimation" "help phdid_postestimation"}{...}
{viewerjumpto "Syntax" "phdid##syntax"}{...}
{viewerjumpto "Description" "phdid##description"}{...}
{viewerjumpto "The problem" "phdid##problem"}{...}
{viewerjumpto "Subcommands" "phdid##subcommands"}{...}
{viewerjumpto "Remarks" "phdid##remarks"}{...}
{viewerjumpto "Examples" "phdid##examples"}{...}
{viewerjumpto "Stored results" "phdid##results"}{...}
{viewerjumpto "References" "phdid##references"}{...}
{viewerjumpto "Authors" "phdid##authors"}{...}
{title:Title}

{phang}
{bf:phdid} {hline 2} Partial homogeneity in staggered difference-in-differences


{marker syntax}{...}
{title:Syntax}

{pstd}Assemble the first-stage estimates{p_end}
{p 8 15 2}
{cmd:phdid setup} [{cmd:,} {it:setup_options}]

{pstd}Ask whether there is heterogeneity to recover{p_end}
{p 8 15 2}
{cmd:phdid test} [{cmd:,} {it:test_options}]

{pstd}Estimate{p_end}
{p 8 15 2}
{cmd:phdid l0} [{cmd:,} {it:l0_options}]

{p 8 15 2}
{cmd:phdid bayes} [{cmd:,} {it:bayes_options}]

{p 8 15 2}
{cmd:phdid fit}{cmd:,} {cmdab:part:ition(}{it:numlist}{cmd:)} [{it:fit_options}]

{pstd}After estimation{p_end}
{p 8 15 2}
{cmd:phdid aggregate} [{it:estimand}] [{cmd:,} {it:options}]

{p 8 15 2}
{cmd:phdid partition} [{cmd:,} {cmd:loss(}{cmdab:vi}|{cmdab:binder}{cmd:)} {cmd:post}]

{p 8 15 2}
{cmd:phdid rhat} [{cmd:,} {cmdab:maxl:ag(}{it:#}{cmd:)}]

{p 8 15 2}
{cmd:phdid sensitivity} [{cmd:,} {it:path_options}]

{p 8 15 2}
{cmd:phdid covcheck} [{cmd:,} {it:options}]

{p 8 15 2}
{cmd:phdid enumerate} [{cmd:,} {it:options}]

{p 8 15 2}
{cmd:phdid graph} [{it:kind}] [{cmd:,} {it:options}]

{pstd}Also{p_end}
{p 8 15 2}
{cmd:phdid describe} [{cmd:,} {cmd:nocells}]

{pstd}
{it:estimand} is {opt overall} (the default), {opt dynamic}, {opt group} or
{opt calendar}.  {it:kind} is {opt coclust}, {opt groups} or {opt bic}.

{pstd}
Each subcommand has its own help file; see
{helpb phdid_setup:phdid setup}, {helpb phdid_test:phdid test},
{helpb phdid_l0:phdid l0}, {helpb phdid_bayes:phdid bayes} and
{helpb phdid_postestimation:phdid postestimation}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:phdid} implements the two estimators, the specification tests and the
diagnostics of Arora and Wagle (2026) for staggered difference-in-differences
designs in which the cohort-time average treatment effects are partially
homogeneous -- equal within groups of cells, but not across them.

{pstd}
In practice you run {helpb phdid_setup:phdid setup} once to read in a first
stage, {helpb phdid_test:phdid test} to ask whether the data can support a
grouping at all, and then {helpb phdid_l0:phdid l0} for a single grouping to
report or {helpb phdid_bayes:phdid bayes} for inference that does not assume
the grouping is known.  Everything else is post-estimation.

{pstd}
The package does not estimate treatment effects from micro data.  It takes a
vector of cohort-time effects and their joint covariance from a
heterogeneity-robust first stage -- {helpb csdid}, {helpb jwdid},
{helpb did_imputation}, or any command that posts {cmd:e(b)} and {cmd:e(V)} --
and asks which of those effects the data say are the same.  Lemma 1 of the
paper is what makes this legitimate: everything the estimators need enters
only through {cmd:R'}{bf:V}{cmd:^-1}{bf:b} and {cmd:R'}{bf:V}{cmd:^-1}{cmd:R},
so the first stage is a sufficient statistic and any two routes to the same
{cmd:(b, V)} give the same answer.


{marker problem}{...}
{title:The problem}

{pstd}
With {it:G} treated cohorts and {it:T} periods there are {it:K} cohort-time
cells.  Estimating all {it:K} effects separately is unbiased but inefficient
when some are in fact equal.  Pooling them into a single two-way fixed-effects
coefficient is efficient but biased whenever the heterogeneity is genuine.
Neither corner is right in general, and the choice between them is a
{it:partition-selection} problem on the cells.

{pstd}
Under the true partition the grouped estimator is unbiased and strictly more
efficient than the flexible one -- in a balanced design, exactly {it:|C|}
times more efficient for every cell in a group of size {it:|C|}.  The
partition is of course unknown, and the package recovers it two ways:

{synoptset 30 tabbed}{...}
{pstd}{bf:phdid l0}{p_end}
{synoptline}
{synopt:Returns}one partition{p_end}
{synopt:Selects by}BIC on the agglomeration path, or a fixed lambda{p_end}
{synopt:Intervals}condition on the selected partition{p_end}
{synopt:Coverage, partition uncertain}0.57-0.62{p_end}
{synopt:Coverage, effects separated}0.95{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 30 tabbed}{...}
{pstd}{bf:phdid bayes}{p_end}
{synoptline}
{synopt:Returns}a posterior over partitions{p_end}
{synopt:Selects by}Dirichlet process prior, collapsed Gibbs{p_end}
{synopt:Intervals}marginalise over the partition{p_end}
{synopt:Coverage, partition uncertain}0.79-0.81{p_end}
{synopt:Coverage, effects separated}0.94{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
The two are not rivals.  {cmd:phdid l0} is the fixed-variance MAP of the same
Bayesian model, differing only in the partition prior.  Use {cmd:phdid l0}
when you want one interpretable grouping to report, and {cmd:phdid bayes} for
inference.


{marker subcommands}{...}
{title:Subcommands}

{pstd}
These are listed in the order you would use them: read in a first stage, ask
whether there is heterogeneity to recover, estimate, then report.

{synoptset 24 tabbed}{...}
{synopthdr:subcommand}
{synoptline}
{syntab:Input}
{synopt:{helpb phdid_setup:setup}}read the first stage and keep the
post-treatment cells{p_end}
{synopt:{opt describe}}report what is in memory{p_end}

{syntab:Specification}
{synopt:{helpb phdid_test:test}}is there heterogeneity to recover?{p_end}

{syntab:Estimation}
{synopt:{helpb phdid_l0:l0}}the l0-penalised estimator, Appendix B{p_end}
{synopt:{helpb phdid_bayes:bayes}}the Dirichlet process estimator{p_end}
{synopt:{opt fit}}restricted GLS under a partition you supply{p_end}

{syntab:After estimation}
{synopt:{opt aggregate}}overall, event-study, by-cohort, by-period{p_end}
{synopt:{opt partition}}a point partition from the posterior{p_end}
{synopt:{opt rhat}}Gelman-Rubin statistics and effective sample sizes{p_end}
{synopt:{opt sensitivity}}regularisation paths over alpha or lambda{p_end}
{synopt:{opt covcheck}}exact against diagonal covariance, Remark 1{p_end}
{synopt:{opt enumerate}}the exact partition posterior, Appendix E{p_end}
{synopt:{opt graph}}co-clustering, grouped effects, the selection path{p_end}
{synoptline}
{p2colreset}{...}


{marker remarks}{...}
{title:Remarks}

{dlgtab:The covariance is not diagonal}

{pstd}
Cohort-time cells share units and share the unit and time fixed effects, so
their first-stage estimates are correlated.  Using the exact joint covariance
is what makes the reported intervals honest.  Treating the cells as
independent understates the posterior variance of aggregates by about 3.5
times in the paper's design, and can misstate individual co-clustering
probabilities by up to 0.39.  {cmd:phdid setup} warns when it is handed a
diagonal covariance, and {cmd:phdid covcheck} quantifies the difference on
your own design.  This is Remark 1 of the paper.

{dlgtab:A partition is only worth reporting when it can be recovered}

{pstd}
Both estimators will partition pure noise if asked to.  Below the separation
threshold the groups they return are quantiles of sampling error, and the
output does not look any different.  In the paper's second application the l0
estimator splits 138 event-study effects into five neat bands while a
randomization test cannot reject a single common effect.  Run
{helpb phdid_test:phdid test} first; it states a reading in words rather than
leaving the output to speak for itself.

{dlgtab:Selection and inference}

{pstd}
The deviance reported by {cmd:phdid l0} is not a valid test of the selected
partition, because the partition was chosen to make it small.  The same
statistic from {cmd:phdid fit}, on a partition fixed in advance, is a valid
test.  More generally, the intervals from {cmd:phdid l0} condition on the
selected partition being correct and under-cover when it is uncertain; that
is what {cmd:phdid bayes} is for.


{marker examples}{...}
{title:Examples}

{pstd}
What follows is the paper's first application -- the effect of state
minimum-wage increases on county-level teen employment, using the panel of
Callaway and Sant'Anna (2021) -- worked through in the order the paper
presents it.  Each block is one step of that argument, and the numbers it
prints are the numbers reported in the paper, so running it end to end also
confirms the installation is sound.  The same sequence, with fuller comments,
is {bf:examples/mpdta.do} in the package source.

{pstd}
{bf:The first stage.}  {cmd:phdid} does not touch the micro data.
{cmd:csdid} estimates the {it:K} cohort-time effects and their joint
covariance, and {cmd:phdid setup} reads those out of {cmd:e(b)} and
{cmd:e(V)} and keeps the post-treatment cells.  Any first stage that posts
the two would serve here.{p_end}

{phang2}{cmd:. webuse set https://friosavila.github.io/playingwithstata/drdid/}{p_end}
{phang2}{cmd:. webuse mpdta, clear}{p_end}
{phang2}{cmd:. csdid lemp, ivar(countyreal) time(year) gvar(first_treat) method(dripw) notyet}{p_end}
{phang2}{cmd:. phdid setup}{p_end}

{pstd}
{bf:Is there anything to recover?}  Both estimators will partition pure noise
if asked to, so this question comes before any estimation.  The test compares
the flexible fit against a single common effect and reports how much of the
spread across cells is signal rather than sampling error.{p_end}

{phang2}{cmd:. phdid test}{p_end}

{pstd}
{bf:One grouping to report.}  {cmd:phdid l0} walks the agglomeration path and
selects a partition by BIC.  {cmd:phdid aggregate} collapses the grouped cell
effects into the overall ATT, and {cmd:phdid graph groups} draws the grouped
estimate against the flexible one for each cell, which is where the gain in
precision becomes visible.{p_end}

{phang2}{cmd:. phdid l0}{p_end}
{phang2}{cmd:. phdid aggregate}{p_end}
{phang2}{cmd:. phdid graph groups}{p_end}

{pstd}
{bf:Inference that does not assume the grouping is known.}  The intervals
above condition on the selected partition being the right one.
{cmd:phdid bayes} samples partitions instead, so its intervals carry the
uncertainty about the grouping itself.  {cmd:phdid rhat} checks that the four chains have
mixed before any of it is believed; {cmd:phdid aggregate dynamic} reports the
event-study profile, averaging cells by time since treatment; and
{cmd:phdid graph coclust} shows, for every pair of cells, the posterior
probability that the two belong to the same group.{p_end}

{phang2}{cmd:. set seed 20260301}{p_end}
{phang2}{cmd:. phdid bayes, alpha(1) iters(20000) burnin(2000) chains(4) init(dispersed)}{p_end}
{phang2}{cmd:. phdid rhat}{p_end}
{phang2}{cmd:. phdid aggregate dynamic}{p_end}
{phang2}{cmd:. phdid graph coclust}{p_end}

{pstd}
{bf:How much of that is the prior talking?}  The concentration parameter
{cmd:alpha} governs how readily the prior opens a new group, and it is not
estimated from the data.  This traces every cell's estimate across a range of
it; a conclusion that holds across the whole range is not an artefact of the
prior.{p_end}

{phang2}{cmd:. phdid sensitivity, alpha(0.25 0.5 1 2 4 8) type(cells)}{p_end}

{pstd}
{bf:The exact posterior, for comparison.}  With {it:K} this small all 877
partitions of the cells can be enumerated, so the posterior can be computed
without sampling at all.  This is the benchmark the sampler is certified
against, and on a design of this size it is the answer to prefer.{p_end}

{phang2}{cmd:. phdid enumerate, alpha(1)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:phdid l0}, {cmd:phdid bayes} and {cmd:phdid fit} are {cmd:eclass} and
post {cmd:e(b)} and {cmd:e(V)} over the {it:K} cohort-time effects, so
{helpb lincom}, {helpb nlcom}, {helpb test} and {helpb estimates} work on
them in the usual way.  The other subcommands are {cmd:rclass}.  See each
subcommand's help file for the full list.

{pstd}
For {cmd:phdid bayes}, {cmd:e(V)} is the posterior covariance of the cell
effects across draws, which already includes uncertainty about the partition.
For {cmd:phdid l0} and {cmd:phdid fit} it is the restricted GLS covariance,
which does not.


{marker references}{...}
{title:References}

{phang}
Arora, P., and R. Wagle. 2026. A Bayesian Approach to Partial Homogeneity in
Staggered Difference-in-Difference. Ashoka University Economics Discussion
Paper 166.
{browse "https://www.ashoka.edu.in/research/a-bayesian-approach-to-partial-homogeneity-in-staggered-difference-in-difference/":link}

{phang}
Arora, P., and R. Wagle. 2026. A Bayesian Approach to Partial Homogeneity in
Staggered Difference-in-Difference (March 1, 2026). SSRN Working Paper
7207083. {browse "https://doi.org/10.2139/ssrn.7207083":doi:10.2139/ssrn.7207083}

{phang}
Callaway, B., and P. H. C. Sant'Anna. 2021. Difference-in-Differences with
Multiple Time Periods. {it:Journal of Econometrics} 225(2): 200-230.

{phang}
Neal, R. M. 2000. Markov Chain Sampling Methods for Dirichlet Process Mixture
Models. {it:Journal of Computational and Graphical Statistics} 9(2): 249-265.

{phang}
Wade, S., and Z. Ghahramani. 2018. Bayesian Cluster Analysis: Point Estimation
and Credible Balls. {it:Bayesian Analysis} 13(2): 559-626.

{phang}
Wooldridge, J. M. 2025. Two-Way Fixed Effects, the Two-Way Mundlak Regression,
and Difference-in-Differences Estimators. {it:Empirical Economics} 69:
2545-2587.


{marker citation}{...}
{title:Citation}

{pstd}
If you use {cmd:phdid}, please cite {bf:both} the Ashoka University discussion
paper and the SSRN working paper listed above, and the software itself.

{pstd}
An R package of the same name, with the same estimators, is on CRAN; the Stata
implementation is certified against it to machine precision on the paper's
first application.  See {bf:cert/phdid_cert.do} in the package source.


{marker authors}{...}
{title:Authors}

{pstd}
The methods implemented here are joint work by Parush Arora and Rohan Wagle,
Department of Economics, Ashoka University.

{pstd}
Rohan Wagle (maintainer){break}
Department of Economics, Ashoka University{break}
rohan.wagle_ug2023@ashoka.edu.in

{pstd}
Parush Arora{break}
Department of Economics, Ashoka University{break}
parush.arora@ashoka.edu.in

{pstd}
Distributed under the MIT license, Copyright (c) 2026 Rohan Wagle and Parush
Arora.
