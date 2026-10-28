{smcl}
{* *! version 0.1.0  26oct2026}{...}
{vieweralsosee "phdid" "help phdid"}{...}
{vieweralsosee "phdid test" "help phdid_test"}{...}
{vieweralsosee "phdid l0" "help phdid_l0"}{...}
{vieweralsosee "phdid bayes" "help phdid_bayes"}{...}
{viewerjumpto "Syntax" "phdid_setup##syntax"}{...}
{viewerjumpto "Description" "phdid_setup##description"}{...}
{viewerjumpto "Options" "phdid_setup##options"}{...}
{viewerjumpto "Remarks" "phdid_setup##remarks"}{...}
{viewerjumpto "Examples" "phdid_setup##examples"}{...}
{viewerjumpto "Stored results" "phdid_setup##results"}{...}
{title:Title}

{phang}
{bf:phdid setup} {hline 2} Assemble the first-stage input


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:phdid setup} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt b(name)}}matrix of first-stage effects; default is {cmd:e(b)}{p_end}
{synopt:{opt v(name)}}their joint covariance; default is {cmd:e(V)}{p_end}
{synopt:{opt cells(string)}}cohort-time names, one per column of {opt b()}{p_end}
{synopt:{opt w:eights(name)}}aggregation weights, one per cell{p_end}
{synopt:{opt nu:nits(#)}}number of units in the first stage{p_end}
{synopt:{opt nowarn}}suppress the diagonal-covariance warning{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{cmd:phdid setup} is {cmd:rclass}.  It must be run before any other
{cmd:phdid} subcommand, and the data it stores persist until it is run again
or Stata's Mata memory is cleared.


{marker description}{...}
{title:Description}

{pstd}
{cmd:phdid setup} reads a vector of cohort-time average treatment effects and
their joint covariance, drops the pre-treatment cells, recovers the cohort and
period of each remaining cell, picks up aggregation weights where the first
stage posts them, and stores the result for the estimation subcommands.

{pstd}
With no options it reads {cmd:e(b)} and {cmd:e(V)} from whatever estimation
command was run last.  This is the usual route: run {helpb csdid},
{helpb jwdid}, {helpb did_imputation} or any other heterogeneity-robust
estimator, then run {cmd:phdid setup}.

{pstd}
Lemma 1 of Arora and Wagle (2026) is what licenses this.  Everything the
partial-homogeneity estimators require enters the likelihood only through
{cmd:R'}{bf:V}{cmd:^-1}{bf:b} and {cmd:R'}{bf:V}{cmd:^-1}{cmd:R}, where
{cmd:R} is the partition indicator.  The first stage is therefore a sufficient
statistic, and two first stages that produce the same {cmd:(b, V)} produce the
same partition posterior, whatever micro data lie behind them.


{marker options}{...}
{title:Options}

{phang}
{opt b(name)} and {opt v(name)} supply the first-stage effects and their
covariance directly, instead of reading {cmd:e(b)} and {cmd:e(V)}.  Both must
be given together.  {opt b()} may be a row or column vector; {opt v()} must be
symmetric and positive definite with as many rows as {opt b()} has elements.
This is the route for effects computed elsewhere, read from a file, or
produced by a command that does not post them in the usual place.

{phang}
{opt cells(string)} gives the cohort-time name of every column of {opt b()},
as a space-separated list such as {cmd:cells(2004:2004 2004:2005 2006:2006)}.
Use it when the coefficient names cannot be parsed.  The names supplied are
themselves parsed, so {cmd:2004:2005} still yields cohort 2004 and period
2005, and cells with period before cohort are still dropped.

{phang}
{opt weights(name)} supplies the weights used to form aggregates, one per
post-treatment cell, as a row or column vector.  They are rescaled to sum to
one.  When this is not given, {cmd:phdid setup} looks for cohort shares posted
by the first stage -- {cmd:csdid} posts them in a {cmd:wgt} equation -- and
falls back to equal weights if there are none.  The weights matter for every
aggregate but for nothing else: the partition posterior does not depend on
them.

{phang}
{opt nunits(#)} records the number of units behind the first stage.  It is
reported by {cmd:phdid describe} and is otherwise unused; it is read from
{cmd:e(N)} when available.

{phang}
{opt nowarn} suppresses the warning issued when the supplied covariance is
diagonal.  Use it only when the cells really are independent, which in a
staggered design they are not.


{marker remarks}{...}
{title:Remarks}

{dlgtab:Which cells are kept}

{pstd}
Only post-treatment cells, those with period at or after the cohort's
treatment date, enter the partition problem.  Pre-treatment cells are the
material for the placebo gauge in {helpb phdid_test:phdid test}, which takes
them as an argument rather than from the stored data.

{dlgtab:Coefficient names}

{pstd}
{cmd:phdid setup} recognises the conventions used by the common first-stage
commands: {cmd:csdid}'s equation names ({cmd:g2004}) with the period as the
final run of digits in the coefficient name; a bare numeric equation name
carrying the cohort, which is how a hand-built matrix with column names
{cmd:2004:2005} arrives, Stata reading the colon as the equation separator;
names of the form {cmd:2004:2005} or {cmd:2004_2005}; and factor-variable
interaction names.  When none matches, every column is kept and the cohort and
period are left unknown, which disables only {cmd:phdid aggregate dynamic},
{cmd:group} and {cmd:calendar}.  Supply {opt cells()} in that case.

{dlgtab:The diagonal warning}

{pstd}
A diagonal covariance is almost always a sign that the cross-cell terms were
discarded somewhere upstream, not that they are zero.  Cohort-time cells share
units and share the unit and time fixed effects.  Running on a diagonal
covariance will produce output, and the output will understate the uncertainty
in every aggregate.  See {cmd:phdid covcheck}.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. csdid lemp, ivar(countyreal) time(year) gvar(first_treat) method(dripw) notyet}{p_end}
{phang2}{cmd:. phdid setup}{p_end}
{phang2}{cmd:. phdid describe}{p_end}

{pstd}From matrices computed elsewhere:{p_end}

{phang2}{cmd:. matrix colnames TAU = 2004 2005 2006 2007 2006 2007 2007}{p_end}
{phang2}{cmd:. matrix coleq    TAU = 2004 2004 2004 2004 2006 2006 2007}{p_end}
{phang2}{cmd:. phdid setup, b(TAU) v(SIGMA)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:phdid setup} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(K)}}number of post-treatment cells{p_end}
{synopt:{cmd:r(isdiag)}}1 if the covariance supplied was diagonal{p_end}
{synopt:{cmd:r(nunits)}}units in the first stage, if known{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(cells)}}cell effects, standard errors and weights{p_end}
{p2colreset}{...}
