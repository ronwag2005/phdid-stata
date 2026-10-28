*! version 0.1.0  phdid describe -- what is in memory
*! Rohan Wagle and Parush Arora
program define phdid_desc, rclass
    version 17.0
    syntax [, noCELLS ]

    mata: phdid_describe()

    tempname cm
    matrix `cm' = r(cells)
    local K      = r(K)
    local isdiag = r(isdiag)
    local nunits = r(nunits)
    local snr    = r(snr_desc)

    di ""
    di as txt "{hline 72}"
    di as txt "phdid: first-stage effects in memory"
    di as txt "{hline 72}"
    di as txt "  Cohort-time cells (K) : " as res `K'
    if (`nunits' < .) {
        di as txt "  Units in first stage  : " as res %12.0fc `nunits'
    }
    di as txt "  Covariance supplied   : " as res ///
        cond(`isdiag', "diagonal", "full (exact cross-cell terms)")
    di as txt "  sd(tau)/median(se)    : " as res %6.2f `snr' ///
        as txt "   (descriptive only; see {bf:phdid test})"
    di as txt "{hline 72}"
    if ("`cells'" == "") {
        matlist `cm', format(%9.5f) twidth(14) ///
            title("  Cell effects, standard errors and aggregation weights")
    }
    if (`isdiag') {
        di ""
        di as err "warning: the covariance in memory is diagonal."
        di as txt "{p 4 4 2}Cohort-time cells share units and share the " ///
            "unit and time fixed effects, so their first-stage estimates " ///
            "are correlated.  Discarding those terms understates the " ///
            "posterior variance of aggregates and misstates co-clustering " ///
            "probabilities; see {bf:phdid covcheck} and Remark 1 of the " ///
            "paper.{p_end}"
    }

    return scalar K = `K'
    return scalar isdiag = `isdiag'
    return matrix cells = `cm'
end
