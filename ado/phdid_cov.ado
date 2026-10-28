*! version 0.1.0  phdid covcheck -- exact vs diagonal covariance
*! Rohan Wagle and Parush Arora
program define phdid_cov, rclass
    version 17.0
    syntax [, ALPHA(real 1) MU0(string) S0SQ(string) ITERs(integer 5000) ///
              BURNin(integer 1000) THIN(integer 1) CHAINs(integer 1) ///
              INIT(string) SEED(string) ]

    if ("`seed'" != "") set seed `seed'
    if ("`init'" == "") local init singletons
    local mu0v = cond("`mu0'" == "", ".", "`mu0'")
    local s0v  = cond("`s0sq'" == "", ".", "`s0sq'")

    di as txt "running the sampler twice (exact, then diagonal) ..."
    mata: phdid_covcheck_run(`alpha', `mu0v', `s0v', `iters', `burnin', ///
                             `thin', `chains', "`init'")

    tempname cm
    matrix `cm' = r(cells)
    local mx = r(cc_maxdiff)
    local vr = r(var_ratio)
    local isdiag = r(isdiag)

    di ""
    di as txt "{hline 72}"
    di as txt "Exact versus diagonal handling of the first-stage covariance"
    di as txt "{hline 72}"
    di as txt "                                   " _col(38) "exact" _col(52) "diagonal"
    di as txt "  Posterior E[# groups]            " _col(36) as res %8.4f r(m_exact) _col(50) %8.4f r(m_diag)
    di as txt "  Overall effect                   " _col(36) as res %8.4f r(att_exact) _col(50) %8.4f r(att_diag)
    di as txt "  Posterior sd of overall effect   " _col(36) as res %8.5f r(sd_exact) _col(50) %8.5f r(sd_diag)
    di as txt "{hline 72}"
    di as txt "  Variance ratio (exact/diagonal)  : " as res %8.2f `vr'
    di as txt "  Largest co-clustering discrepancy: " as res %8.3f `mx'
    di as txt "{hline 72}"
    if (`isdiag') {
        di as txt "{p 2 2 2}The covariance in memory is itself diagonal, so " ///
            "the two runs differ only by Monte Carlo error.  To use this " ///
            "check as intended, supply the exact joint covariance of the " ///
            "first-stage effects.{p_end}"
    }
    else {
        di as txt "{p 2 2 2}The diagonal shortcut drops the cross-cell " ///
            "terms in the cluster-assignment moves.  A variance ratio " ///
            "above one means it understates the uncertainty in the " ///
            "aggregate by that factor, and the discrepancy above is how " ///
            "far a single co-clustering probability moves.  This is " ///
            "Remark 1 of the paper; the package default is exact.{p_end}"
    }

    return scalar var_ratio = `vr'
    return scalar cc_maxdiff = `mx'
    return matrix cells = `cm'
end
