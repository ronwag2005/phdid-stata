*! version 0.1.0  phdid test -- is there heterogeneity to recover?
*! Rohan Wagle and Parush Arora
program define phdid_test, rclass
    version 17.0
    syntax [, PRE(name) NSIM(integer 10000) noCENTER SEED(string) ]

    if ("`seed'" != "") set seed `seed'
    local ctr = cond("`center'" == "", 1, 0)

    tempname prem
    if ("`pre'" != "") {
        matrix `prem' = `pre'
        mata: phdid_test_run(st_matrix("`prem'")', `nsim', `ctr')
    }
    else {
        mata: phdid_test_run(J(0, 1, .), `nsim', `ctr')
    }

    di ""
    di as txt "{hline 72}"
    di as txt "phdid: is there heterogeneity to recover?"
    di as txt "{hline 72}"
    di as txt "1. Common-effect test (pooled GLS deviance)" _col(50) ///
        as txt "[package addition]"
    di as txt "     chi2(" as res r(df) as txt ") = " as res %8.2f r(chi2) ///
        as txt ",  Prob > chi2 = " as res %8.4f r(p)
    di ""
    di as txt "2. Dispersion decomposition" _col(50) as txt "[paper 5.2.3, adapted]"
    di as txt "     observed cross-cell sd        : " as res %9.5f sqrt(r(disp_obs))
    di as txt "     expected under a common effect: " as res %9.5f sqrt(r(disp_exp))
    di as txt "     heterogeneity sd (excess)     : " as res %9.5f sqrt(r(het_var))
    di as txt "     heterogeneity share           : " as res %8.0f 100*r(het_share) as txt "%"
    di as txt "     signal-to-noise               : " as res %9.2f r(snr)
    di as txt "     sd(tau)/median(se)            : " as res %9.2f r(snr_desc) ///
        as txt "   (descriptive)"

    if (r(has_placebo)) {
        di ""
        di as txt "3. Placebo gauge and randomization test" _col(50) ///
            as txt "[paper 5.2.3]"
        di as txt "     s_pre  (noise)                : " as res %9.5f r(s_pre)
        di as txt "     s_post (observed spread)      : " as res %9.5f r(s_post)
        di as txt "     signal-to-noise               : " as res %9.2f r(pl_snr)
        di as txt "     heterogeneity share           : " as res %8.0f 100*r(pl_het_share) as txt "%"
        di as txt "     randomization p (one-sided)   : " as res %9.3f r(pl_p) ///
            as txt "  (" as res r(nsim) as txt " draws)"
    }

    di as txt "{hline 72}"
    phdid_test_verdict
    di ""
    di as txt "{p 2 2 2}Tags mark provenance: [paper 5.2.3] is the published " ///
        "procedure; [package addition] and [adapted] are this package's. " ///
        "The reading above is the package's rule, not the paper's.{p_end}"

    return scalar chi2 = r(chi2)
    return scalar df   = r(df)
    return scalar p    = r(p)
    return scalar het_share = r(het_share)
    return scalar snr  = r(snr)
end

program define phdid_test_verdict
    version 17.0
    local p     = r(p)
    local share = r(het_share)
    local snr   = r(snr)
    local prand = cond(r(has_placebo), r(pl_p), .)

    local reject = (`p' < 0.05)
    if (`prand' < . & `prand' >= 0.10) local reject 0

    if (`reject' == 0 & `share' < 0.10) {
        di as txt "Reading: no recoverable heterogeneity."
        di as txt "{p 2 2 2}A common effect cannot be rejected and the " ///
            "estimated heterogeneity share is near zero, so pooling is an " ///
            "adequate specification for the aggregate. Groups returned by " ///
            "{bf:phdid l0} or {bf:phdid bayes} here would be quantiles of " ///
            "noise; do not read them as effect clusters. This is a failure " ///
            "to reject, not proof of a literally common effect.{p_end}"
        exit
    }
    if (`reject' == 1 & `snr' >= 2) {
        di as txt "Reading: recoverable heterogeneity."
        di as txt "{p 2 2 2}The spread across cells clearly exceeds what " ///
            "sampling noise would produce, so there is structure for the " ///
            "partition model to find. Proceed to {bf:phdid l0} for a point " ///
            "partition and {bf:phdid bayes} for inference, and check the " ///
            "co-clustering matrix to see which groupings are firm.{p_end}"
        exit
    }
    di as txt "Reading: borderline, the low-separation regime."
    di as txt "{p 2 2 2}There is some evidence of heterogeneity, but the " ///
        "effects are not clearly separated. This is the regime the paper " ///
        "identifies as the boundary of reliable recovery, where selection " ///
        "noise can erase or even reverse the precision gain and the two " ///
        "estimators may disagree. Treat any single partition with caution.{p_end}"
end
