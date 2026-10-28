*! version 0.1.0  phdid rhat -- sampler convergence
*! Rohan Wagle and Parush Arora
program define phdid_rhat, rclass
    version 17.0
    syntax [, MAXLag(integer 200) ]

    if ("`e(cmd)'" != "phdid_bayes") {
        di as err "no phdid bayes results found"
        exit 301
    }

    mata: phdid_rhat_run(`maxlag')

    tempname tbl
    matrix `tbl' = r(table)
    local nch  = r(nchains)
    local rmax = r(rhat_max)
    local emin = r(ess_min)

    di ""
    di as txt "{hline 72}"
    di as txt "Sampler diagnostics" _col(40) as txt "`nch' chain(s), split in half"
    di as txt "{hline 72}"
    matlist `tbl', format(%9.3f) twidth(14)
    di as txt "{hline 72}"
    if (`nch' < 2) {
        di as txt "{p 2 2 2}With a single chain these are split-half " ///
            "statistics only.  Re-run with {bf:chains(4) init(dispersed)} " ///
            "for a diagnostic that can detect a chain stuck in one corner " ///
            "of partition space.{p_end}"
    }
    else if (`rmax' > 1.05 | `emin' < 400) {
        di as txt "Reading: the sampler has not clearly converged."
        di as txt "{p 2 2 2}Rhat above 1.05 or an effective sample below " ///
            "about 400 means the reported intervals are not yet reliable. " ///
            "Increase {bf:iters()}, and check whether the chains disagree " ///
            "about the number of groups rather than only about the cell " ///
            "effects.{p_end}"
    }
    else {
        di as txt "Reading: converged."
        di as txt "{p 2 2 2}The chains agree on both the cell effects and " ///
            "the number of groups, so the posterior summaries can be " ///
            "read as reported.{p_end}"
    }

    return matrix table = `tbl'
    return scalar rhat_max = `rmax'
    return scalar ess_min  = `emin'
    return scalar nchains  = `nch'
end
