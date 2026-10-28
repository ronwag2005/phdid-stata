*! version 0.1.0  phdid: the paper's first application
*! Rohan Wagle and Parush Arora
*
*  Reproduces the first application of Arora and Wagle (2026) on the
*  county-level minimum-wage panel of Callaway and Sant'Anna (2021),
*  exercising every subcommand in the order you would actually use them.
*
*      do examples/mpdta.do
*
*  Requires csdid for the first stage:
*      ssc install csdid
*      ssc install drdid

version 17.0
clear all
set more off
set linesize 100

* ------------------------------------------------------------------
* 0. the first stage
*
*    Any heterogeneity-robust estimator will do.  phdid reads e(b) and
*    e(V), which Lemma 1 says is all it needs.
* ------------------------------------------------------------------
use "https://friosavila.github.io/playingwithstata/drdid/mpdta.dta", clear
csdid lemp, ivar(countyreal) time(year) gvar(first_treat) method(dripw) notyet

phdid setup
phdid describe

* ------------------------------------------------------------------
* 1. is there heterogeneity to recover?
*
*    Run this before either estimator.  Both will partition pure noise
*    if asked to, and the output looks no different when they do.
* ------------------------------------------------------------------
phdid test

* ------------------------------------------------------------------
* 2. one partition to report
* ------------------------------------------------------------------
phdid l0

phdid aggregate
phdid aggregate dynamic
phdid aggregate group

phdid graph groups, name(l0groups, replace)
phdid graph bic,    name(l0path,   replace)

* A partition fixed in advance, rather than selected: here the deviance
* is an honest test of the restriction.
phdid fit, partition(1 1 2 2 3 3 3)

* ------------------------------------------------------------------
* 3. inference that does not assume the partition is known
* ------------------------------------------------------------------
set seed 20260301
phdid bayes, alpha(1) iters(20000) burnin(2000) chains(4) init(dispersed)

phdid rhat
phdid aggregate
phdid aggregate dynamic

phdid partition, loss(vi)
phdid graph coclust, name(coclust, replace)

* ------------------------------------------------------------------
* 4. how much of this is the prior talking?
*
*    The aggregate is the quantity least sensitive to over-pooling, so
*    type(cells) is the view that shows where the regularisation bites.
* ------------------------------------------------------------------
phdid sensitivity, alpha(0.25 0.5 1 2 4 8) type(overall) ///
    iters(6000) burnin(1000) name(alphapath, replace)

phdid sensitivity, by(m) type(cells) name(mpath, replace)

* ------------------------------------------------------------------
* 5. does the exact covariance matter here?  (Remark 1)
* ------------------------------------------------------------------
set seed 7
phdid covcheck, alpha(1) iters(6000) burnin(1000)

* ------------------------------------------------------------------
* 6. with only seven cells the posterior can be computed exactly
*
*    877 partitions.  This is the benchmark the sampler above should
*    have reproduced; compare r(m_mean) with e(m_mean) from step 3.
* ------------------------------------------------------------------
phdid enumerate, alpha(1)
