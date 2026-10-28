{smcl}
{* *! version 0.1.0  26oct2026}{...}
{vieweralsosee "phdid" "help phdid"}{...}
{vieweralsosee "phdid l0" "help phdid_l0"}{...}
{vieweralsosee "phdid bayes" "help phdid_bayes"}{...}
{viewerjumpto "Syntax" "phdid_test##syntax"}{...}
{viewerjumpto "Description" "phdid_test##description"}{...}
{viewerjumpto "Options" "phdid_test##options"}{...}
{viewerjumpto "What is tested" "phdid_test##what"}{...}
{viewerjumpto "Provenance" "phdid_test##provenance"}{...}
{viewerjumpto "Examples" "phdid_test##examples"}{...}
{viewerjumpto "Stored results" "phdid_test##results"}{...}
{title:Title}

{phang}
{bf:phdid test} {hline 2} Is there heterogeneity to recover?


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:phdid test} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt pre(name)}}vector of pre-treatment (placebo) effects{p_end}
{synopt:{opt nsim(#)}}draws for the randomization test; default
{cmd:nsim(10000)}{p_end}
{synopt:{opt nocenter}}do not centre the placebo effects before resampling{p_end}
{synopt:{opt seed(string)}}set the random-number seed first{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:phdid test} asks whether the cohort-time effects are separated enough for
a partition to be recoverable, before either estimator is run.  This is not a
formality.  Both {cmd:phdid l0} and {cmd:phdid bayes} will partition pure noise
if asked to, and their output looks no different when they do.  In the second
application of Arora and Wagle (2026) the l0 estimator splits 138 event-study
effects into five tidy bands while a randomization test cannot reject a single
common effect.

{pstd}
The command reports three things and then states a reading in words.


{marker what}{...}
{title:What is tested}

{dlgtab:1. The common-effect test}

{pstd}
The deviance of the fully pooled restricted GLS fit, referred to a chi-squared
distribution on {it:K}-1 degrees of freedom.  Under a common effect and the
first-stage covariance taken as known, that deviance is the likelihood-ratio
statistic against the flexible fit.  A small p-value says a single effect does
not describe the cells.

{dlgtab:2. The dispersion decomposition}

{pstd}
The observed cross-cell dispersion of the first-stage effects is compared with
what sampling noise alone would produce under a common effect, computed from
the exact joint covariance rather than from the diagonal.  The excess is the
heterogeneity variance; the heterogeneity share is its fraction of the
observed dispersion, and the signal-to-noise ratio is the ratio of the two
standard deviations.  A share near zero means the spread you see is the spread
you would see anyway.

{dlgtab:3. The placebo gauge and randomization test}

{pstd}
Reported only when {opt pre()} supplies pre-treatment estimates.  Their
dispersion {it:s_pre} is a direct measurement of the noise, with no model
behind it.  Comparing it with the post-treatment dispersion {it:s_post} gives a
model-free signal-to-noise ratio, and resampling the centred placebo effects
gives a one-sided p-value for the observed post-treatment spread.  This is the
procedure of section 5.2.3 of the paper, and it is the most credible of the
three because it does not rely on the covariance being right.


{marker options}{...}
{title:Options}

{phang}
{opt pre(name)} supplies a vector of pre-treatment or placebo cohort-time
effects, as a row or column vector.  These are not taken from the stored data,
because {cmd:phdid setup} drops them: pass them explicitly, for example from
the pre-treatment columns of the first stage's {cmd:e(b)}.

{phang}
{opt nsim(#)} sets the number of resamples for the randomization test.  The
default, 10,000, resolves p-values to about 0.01.

{phang}
{opt nocenter} resamples the placebo effects as they are rather than after
removing their mean.  The default centres them, which is the right choice when
the placebos may carry a common pre-trend that is not part of the noise.

{phang}
{opt seed(string)} sets the random-number seed before resampling, so the
reported p-value is reproducible.


{marker provenance}{...}
{title:Provenance}

{pstd}
The output tags each block.  {bf:[paper 5.2.3]} marks the placebo gauge and
randomization test, which are published procedures.  {bf:[adapted]} marks the
dispersion decomposition, which is the paper's idea computed against the exact
covariance rather than the diagonal.  {bf:[package addition]} marks the
common-effect chi-squared test, which follows from the paper's model but is
not reported in it.

{pstd}
The reading printed at the end -- "no recoverable heterogeneity", "recoverable
heterogeneity", or "borderline, the low-separation regime" -- is the
{it:package's} rule, not the paper's, and is stated as such.  It combines the
p-value, the heterogeneity share and the signal-to-noise ratio, and defers to
the randomization test when one is available.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. phdid test}{p_end}

{pstd}With a placebo gauge built from the pre-treatment cells:{p_end}

{phang2}{cmd:. matrix PRE = e(b)[1, "g2004:t_2003_2004".."g2004:t_2003_2005"]}{p_end}
{phang2}{cmd:. phdid test, pre(PRE) nsim(20000) seed(1234)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:phdid test} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(chi2)}}common-effect deviance{p_end}
{synopt:{cmd:r(df)}}its degrees of freedom, {it:K}-1{p_end}
{synopt:{cmd:r(p)}}its p-value{p_end}
{synopt:{cmd:r(het_share)}}heterogeneity share of observed dispersion{p_end}
{synopt:{cmd:r(snr)}}signal-to-noise ratio{p_end}
{p2colreset}{...}
