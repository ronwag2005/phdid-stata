*! version 0.1.0  phdid certification script
*! Rohan Wagle and Parush Arora
*
*  Certifies the Stata implementation against the R package phdid 0.1.0
*  (CRAN), using the first application of Arora and Wagle (2026) as the
*  common input.  The fixtures in cert/fixtures were written by the R
*  package at 17 significant figures and are not recomputed here.
*
*  Everything the paper's estimators do is deterministic given (tau,
*  Sigma) except the Gibbs sampler, and Stata's random number generator
*  is not R's.  So the deterministic parts are certified against R to
*  machine precision, and the sampler is certified against the exact
*  partition posterior from full enumeration -- which both languages
*  compute deterministically, and which is the sampler's own target.
*
*  Usage:   do cert/phdid_cert.do
*  Passes quietly; any failure prints the comparison and exits 9.

version 17.0
clear all
set more off
set linesize 100
* the fixtures carry 17 significant figures; import delimited would
* otherwise store them as floats and throw away nine of them
set type double

local here "`c(pwd)'"
adopath ++ "`here'/ado"
adopath ++ "`here'/src"

* Certify what is shipped.  The compiled library is what an installed
* copy runs from; the .mata sources are only a fallback for a working
* tree where build.do has not been run yet.
capture confirm file "`here'/src/lphdid.mlib"
if (_rc == 0) {
    mata: mata mlib index
    di as txt "certifying the compiled library src/lphdid.mlib"
}
else {
    foreach f in core drivers bayes post {
        qui do "`here'/src/phdid_`f'.mata"
    }
    di as err "src/lphdid.mlib not found; certifying the Mata sources instead"
    di as txt "run {bf:do build.do} and certify again before releasing"
}

* ------------------------------------------------------------------
* assertion helper
* ------------------------------------------------------------------
capture program drop phdid_ck
program define phdid_ck
    version 17.0
    args got want tol name
    if ("`tol'" == "") local tol 1e-8
    local d = abs(`got' - `want')
    local s = max(1, abs(`want'))
    if (`d' / `s' > `tol' | missing(`got')) {
        di as err "FAIL  `name'"
        di as err "      got  " %21.17g `got'
        di as err "      want " %21.17g `want'
        di as err "      diff " %12.3e `d'
        global PHDID_FAIL = $PHDID_FAIL + 1
    }
    else {
        global PHDID_OK = $PHDID_OK + 1
    }
end

global PHDID_FAIL = 0
global PHDID_OK   = 0

* ------------------------------------------------------------------
* fixtures
* ------------------------------------------------------------------
qui import delimited using "`here'/cert/fixtures/mpdta_expected.csv", ///
    clear varnames(1) stringcols(1)
local nfix = _N
forvalues i = 1/`nfix' {
    local k = key[`i']
    local F_`k' = value[`i']
}

qui import delimited using "`here'/cert/fixtures/mpdta_tau.csv", ///
    clear varnames(1) stringcols(1)
local K = _N
matrix TAU = J(1, `K', .)
matrix WGT = J(`K', 1, .)
local cn ""
local ce ""
forvalues i = 1/`K' {
    matrix TAU[1, `i'] = tau[`i']
    matrix WGT[`i', 1] = w[`i']
    local cn `cn' `=t[`i']'
    local ce `ce' `=g[`i']'
}
matrix colnames TAU = `cn'
matrix coleq    TAU = `ce'

* the sigma fixture carries a label column and a label header row, so it
* is read without names and the first row and column are skipped
qui import delimited using "`here'/cert/fixtures/mpdta_sigma.csv", ///
    clear varnames(nonames) rowrange(2)
matrix SIG = J(`K', `K', .)
forvalues i = 1/`K' {
    forvalues j = 1/`K' {
        local vj = `j' + 1
        matrix SIG[`i', `j'] = v`vj'[`i']
    }
}
matrix colnames SIG = `cn'
matrix coleq    SIG = `ce'
matrix rownames SIG = `cn'
matrix roweq    SIG = `ce'

di as txt _n "phdid certification: `nfix' reference values from phdid 0.1.0 (R)"
di as txt "{hline 70}"

* ==================================================================
* 1. setup, and the Lemma 1 claim that the route in does not matter
* ==================================================================
qui phdid setup, b(TAU) v(SIG) weights(WGT) nowarn
phdid_ck `=r(K)' `K' . "setup: K"

qui phdid describe
tempname CELLS
matrix `CELLS' = r(cells)
forvalues i = 1/`K' {
    phdid_ck `=`CELLS'[`i', 1]' `=TAU[1, `i']' 1e-12 "setup: tau[`i']"
    phdid_ck `=`CELLS'[`i', 2]' `=sqrt(SIG[`i', `i'])' 1e-12 "setup: se[`i']"
}
local wsum = 0
forvalues i = 1/`K' {
    local wsum = `wsum' + `CELLS'[`i', 3]
}
phdid_ck `wsum' 1 1e-12 "setup: weights sum to one"

* ==================================================================
* 2. the l0 estimator: path, selection, fit  (paper Appendix B)
* ==================================================================
qui phdid l0, select(bic)
phdid_ck `=e(m)' `F_l0_m' . "l0: selected m"
phdid_ck `=e(att)' `F_l0_att' 1e-8 "l0: overall ATT"
phdid_ck `=e(att_se)' `F_l0_att_se' 1e-8 "l0: overall ATT se"

tempname PATH L0C
matrix `PATH' = e(path)
matrix `L0C'  = e(cells)
forvalues m = 1/`K' {
    phdid_ck `=`PATH'[`m', 2]' `F_l0_dev_`m'' 1e-7 "l0: deviance at m=`m'"
    phdid_ck `=`PATH'[`m', 3]' `F_l0_bic_`m'' 1e-7 "l0: BIC at m=`m'"
    if (`m' < `K') {
        phdid_ck `=`PATH'[`m', 5]' `F_l0_mergecost_`m'' 1e-8 ///
            "l0: merge cost at m=`m'"
        phdid_ck `=`PATH'[`m', 6]' `F_l0_pairprod_`m'' . ///
            "l0: pair product at m=`m'"
    }
}
forvalues i = 1/`K' {
    phdid_ck `=`L0C'[`i', 3]' `F_l0_group_`i'' . "l0: group label of cell `i'"
    phdid_ck `=`L0C'[`i', 4]' `F_l0_tau_`i'' 1e-8 "l0: grouped effect `i'"
    phdid_ck `=`L0C'[`i', 6]' `F_l0_vratio_`i'' 1e-8 "l0: variance ratio `i'"
}

* dynamic aggregation, conditional on the selected partition
qui phdid aggregate dynamic
tempname DYN
matrix `DYN' = r(table)
forvalues e = 1/4 {
    phdid_ck `=`DYN'[`e', 1]' `F_l0_dyn_est_`e'' 1e-8 "l0: event-time est `e'"
    phdid_ck `=`DYN'[`e', 2]' `F_l0_dyn_se_`e''  1e-8 "l0: event-time se `e'"
}

* ==================================================================
* 3. the two corners: fully flexible and fully pooled
* ==================================================================
qui phdid l0, select(m) m(`K')
phdid_ck `=e(deviance)' `F_flex_dev' 1e-6 "flexible: deviance (zero)"
qui phdid aggregate overall
tempname OV
matrix `OV' = r(table)
phdid_ck `=`OV'[1, 1]' `F_flex_att' 1e-8 "flexible: overall ATT"
phdid_ck `=`OV'[1, 2]' `F_flex_att_se' 1e-8 "flexible: overall ATT se"

qui phdid l0, select(m) m(1)
phdid_ck `=e(deviance)' `F_pooled_dev' 1e-7 "pooled: deviance"
tempname PC
matrix `PC' = e(cells)
phdid_ck `=`PC'[1, 4]' `F_pooled_tau' 1e-8 "pooled: common effect"

* a partition supplied by hand, not selected: this is the route that
* lets a grouping from theory be tested against the data
qui phdid fit, partition(1 1 2 2 3 3 3) notable
tempname PH
matrix `PH' = e(cells)
forvalues i = 1/`K' {
    phdid_ck `=`PH'[`i', 4]' `F_phfit_tau_`i'' 1e-10 "ph_fit: effect `i'"
}
phdid_ck `=e(deviance)' `F_phfit_dev' 1e-9 "ph_fit: deviance"
phdid_ck `=e(m)' 3 . "ph_fit: group count"

* ==================================================================
* 4. the specification test  (paper 5.2.3, plus the package's chi2)
* ==================================================================
qui phdid test
phdid_ck `=r(chi2)' `F_chi2' 1e-7 "test: chi2"
phdid_ck `=r(df)' `F_chi2_df' . "test: df"
phdid_ck `=r(p)' `F_chi2_p' 1e-7 "test: p-value"
phdid_ck `=r(het_share)' `F_disp_share' 1e-8 "test: heterogeneity share"
phdid_ck `=r(snr)' `F_disp_snr' 1e-8 "test: signal-to-noise"

* ==================================================================
* 5. exact partition posterior  (paper Appendix E)
*
*    This is the oracle.  Both languages enumerate all 877 partitions
*    of 7 cells deterministically, so it certifies the collapsed
*    marginal likelihood and the CRP prior exactly -- the two pieces
*    the Gibbs sampler is built from.
* ==================================================================
qui phdid enumerate, alpha(1)
phdid_ck `=r(n_partitions)' `F_enum_n' . "enumerate: Bell(7)"
phdid_ck `=r(m_mean)' `F_enum_Em' 1e-8 "enumerate: E[groups]"
phdid_ck `=r(att)' `F_enum_att' 1e-8 "enumerate: posterior mean ATT"
phdid_ck `=r(att_lo)' `F_enum_lo' 1e-8 "enumerate: ATT lower"
phdid_ck `=r(att_hi)' `F_enum_hi' 1e-8 "enumerate: ATT upper"

tempname ECC
matrix `ECC' = r(coclust)
forvalues i = 1/`K' {
    forvalues j = 1/`K' {
        phdid_ck `=`ECC'[`i', `j']' `F_enum_cc_`i'_`j'' 1e-8 ///
            "enumerate: co-clustering (`i',`j')"
    }
}

* ==================================================================
* 6. the Gibbs sampler, against the oracle above
*
*    Stata's RNG is not R's, so the draws cannot match draw for draw.
*    What must hold is that the sampler recovers the exact posterior it
*    targets, to within Monte Carlo error.
* ==================================================================
set seed 20260301
qui phdid bayes, alpha(1) iters(40000) burnin(5000) chains(4) init(dispersed)
phdid_ck `=e(m_mean)' `F_enum_Em' 0.01 "sampler: E[groups] vs oracle"
qui phdid aggregate overall
matrix `OV' = r(table)
phdid_ck `=`OV'[1, 1]' `F_enum_att' 0.02 "sampler: ATT vs oracle"

tempname SCC
matrix `SCC' = e(coclust)
local worst = 0
forvalues i = 1/`K' {
    forvalues j = 1/`K' {
        local d = abs(`SCC'[`i', `j'] - `ECC'[`i', `j'])
        if (`d' > `worst') local worst = `d'
    }
}
phdid_ck `worst' 0 0.02 "sampler: worst co-clustering error vs oracle"

qui phdid rhat
if (r(rhat_max) > 1.02) {
    di as err "FAIL  sampler: max Rhat = " %6.4f r(rhat_max)
    global PHDID_FAIL = $PHDID_FAIL + 1
}
else {
    global PHDID_OK = $PHDID_OK + 1
}

* ==================================================================
* 7. invariants that do not depend on R
* ==================================================================
* relabelling a partition must not change the fit: the group numbers
* are names, so permuting them is a no-op
qui phdid fit, partition(1 1 2 2 3 3 3) notable
local attA = e(att)
local devA = e(deviance)
qui phdid fit, partition(3 3 1 1 2 2 2) notable
phdid_ck `=e(att)' `attA' 1e-14 "invariance: ATT under relabelling"
phdid_ck `=e(deviance)' `devA' 1e-14 "invariance: deviance under relabelling"

* the selected partition cannot fit worse than any coarser one
qui phdid l0, select(bic)
local devsel = e(deviance)
qui phdid fit, partition(1 1 1 1 1 1 1) notable
if (e(deviance) < `devsel' - 1e-10) {
    di as err "FAIL  invariance: pooled deviance below the selected one"
    global PHDID_FAIL = $PHDID_FAIL + 1
}
else {
    global PHDID_OK = $PHDID_OK + 1
}

* aggregation weights are a partition of unity within each term
qui phdid aggregate group
tempname AA
matrix `AA' = r(A)
forvalues i = 1/`=rowsof(`AA')' {
    local rs = 0
    forvalues j = 1/`K' {
        local rs = `rs' + `AA'[`i', `j']
    }
    phdid_ck `rs' 1 1e-12 "aggregate: group row `i' sums to one"
}

* grouping can only shrink the sampling variance of a cell
qui phdid l0, select(bic)
matrix `L0C' = e(cells)
forvalues i = 1/`K' {
    if (`L0C'[`i', 6] > 1 + 1e-10) {
        di as err "FAIL  invariance: cell `i' variance ratio exceeds one"
        global PHDID_FAIL = $PHDID_FAIL + 1
    }
    else {
        global PHDID_OK = $PHDID_OK + 1
    }
}

* the diagonal shortcut really is different on this design (Remark 1)
set seed 7
qui phdid covcheck, alpha(1) iters(4000) burnin(1000)
if (r(cc_maxdiff) < 0.05) {
    di as err "FAIL  covcheck: exact and diagonal agree, expected a gap"
    global PHDID_FAIL = $PHDID_FAIL + 1
}
else {
    global PHDID_OK = $PHDID_OK + 1
}

* ==================================================================
di as txt "{hline 70}"
if ($PHDID_FAIL == 0) {
    di as txt "phdid certification: " as res "$PHDID_OK" as txt " checks, all passed."
    di as txt "{hline 70}"
}
else {
    di as err "phdid certification: $PHDID_FAIL of " ///
        `=$PHDID_OK + $PHDID_FAIL' " checks FAILED."
    exit 9
}
