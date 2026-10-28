*! version 0.1.0  phdid enumerate -- exact partition posterior
*! Rohan Wagle and Parush Arora
program define phdid_enum, rclass
    version 17.0
    syntax [, ALPHA(real 1) MU0(string) S0SQ(string) MAXCells(integer 11) ///
              Level(cilevel) noCOCLUST ]

    local mu0v = cond("`mu0'" == "", ".", "`mu0'")
    local s0v  = cond("`s0sq'" == "", ".", "`s0sq'")

    mata: phdid_enum_run(`alpha', `mu0v', `s0v', `maxcells', `level')

    tempname ccm taum probm
    matrix `ccm'   = r(coclust)
    matrix `taum'  = r(tau)
    matrix `probm' = r(prob)

    di ""
    di as txt "{hline 72}"
    di as txt "Exact partition posterior (enumeration)"
    di as txt "{hline 72}"
    di as txt "  Cells (K)        : " as res r(K)
    di as txt "  Partitions       : " as res %12.0fc r(n_partitions) ///
        as txt "  (enumerated in full)"
    di as txt "  Concentration    : alpha = " as res %6.4g r(alpha)
    di as txt ""
    di as txt "  Posterior E[# groups] : " as res %7.4f r(m_mean)
    di as txt "  Overall effect        : " as res %8.4f r(att) ///
        as txt "  [" as res %7.4f r(att_lo) as txt ", " ///
        as res %7.4f r(att_hi) as txt "]"
    di as txt "{hline 72}"
    if ("`coclust'" == "") {
        di as txt "  Co-clustering probabilities:"
        matlist `ccm', format(%6.2f) twidth(12)
    }
    di ""
    di as txt "{p 2 2 2}This is the deterministic benchmark for " ///
        "{bf:phdid bayes}: the sampler should reproduce these numbers to " ///
        "within Monte Carlo error.{p_end}"

    return scalar n_partitions = r(n_partitions)
    return scalar m_mean = r(m_mean)
    return scalar att    = r(att)
    return scalar att_lo = r(att_lo)
    return scalar att_hi = r(att_hi)
    return matrix coclust = `ccm'
    return matrix tau     = `taum'
end
