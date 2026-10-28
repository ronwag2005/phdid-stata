*! version 0.1.0  phdid bayes -- Dirichlet process partial-homogeneity estimator
*! Rohan Wagle and Parush Arora
program define phdid_bayes, eclass
    version 17.0
    syntax [, ALPHA(real 1) MU0(string) S0SQ(string) ITERs(integer 4000) ///
              BURNin(integer 1000) THIN(integer 1) CHAINS(integer 1)     ///
              MARGinal(string) INIT(string) Level(cilevel) SEED(string)  ///
              noCOCLUST ]

    if ("`marginal'" == "") local marginal "exact"
    if (!inlist("`marginal'", "exact", "diagonal")) {
        di as err "marginal() must be {bf:exact} or {bf:diagonal}"
        exit 198
    }
    if ("`init'" == "") local init "dispersed"
    if (!inlist("`init'", "dispersed", "singletons", "pooled", "random")) {
        di as err "init() must be dispersed, singletons, pooled or random"
        exit 198
    }
    if (`alpha' <= 0) {
        di as err "alpha() must be positive"
        exit 198
    }
    if (`burnin' >= `iters') {
        di as err "burnin() must be smaller than iters()"
        exit 198
    }
    if ("`seed'" != "") set seed `seed'

    * mu0 and s0sq default to data-scaled values chosen in Mata; an empty
    * option is passed through as missing to request that default.
    local mu0v = cond("`mu0'" == "", ".", "`mu0'")
    local s0v  = cond("`s0sq'" == "", ".", "`s0sq'")
    if ("`mu0'" != "") {
        cap confirm number `mu0'
        if (_rc) {
            di as err "mu0() must be a number"
            exit 198
        }
    }
    if ("`s0sq'" != "") {
        cap confirm number `s0sq'
        if (_rc) {
            di as err "s0sq() must be a number"
            exit 198
        }
        if (`s0sq' <= 0) {
            di as err "s0sq() must be positive"
            exit 198
        }
    }

    mata: st_numscalar("r(K)", phdid_K())
    local K = r(K)

    di as txt "  drawing `iters' sweeps x `chains' chain(s) ..."
    mata: phdid_bayes_run(`alpha', `mu0v', `s0v', `iters', `burnin', ///
        `thin', `chains', "`marginal'", "`init'", `level')

    tempname cmat ccmat bmat Vmat
    matrix `cmat'  = r(cells)
    matrix `ccmat' = r(coclust)
    matrix `bmat'  = r(b)
    matrix `Vmat'  = r(V)
    local mmean  = r(m_mean)
    local nd     = r(ndraws)
    local rmu0   = r(mu0)
    local rs0sq  = r(s0sq)
    local priorm = r(prior_m)

    ereturn post `bmat' `Vmat', depname("CATT")
    ereturn local cmd     "phdid_bayes"
    ereturn local cmdline "phdid bayes `0'"
    ereturn local title   "Dirichlet process partial-homogeneity estimator"
    ereturn local marginal "`marginal'"
    ereturn local init     "`init'"
    ereturn matrix cells   = `cmat'
    ereturn matrix coclust = `ccmat'
    ereturn scalar alpha   = `alpha'
    ereturn scalar mu0     = `rmu0'
    ereturn scalar s0sq    = `rs0sq'
    ereturn scalar m_mean  = `mmean'
    ereturn scalar ndraws  = `nd'
    ereturn scalar prior_m = `priorm'
    ereturn scalar chains  = `chains'
    ereturn scalar iters   = `iters'
    ereturn scalar burnin  = `burnin'
    ereturn scalar K       = `K'
    ereturn scalar level   = `level'

    di ""
    di as txt "{hline 72}"
    di as txt "Dirichlet process partial-homogeneity posterior"
    di as txt "{hline 72}"
    di as txt "  Cells (K)        : " as res e(K)
    di as txt "  Concentration    : alpha = " as res %6.4g e(alpha) ///
        as txt "   (prior E[m] = " as res %5.2f e(prior_m) as txt ")"
    di as txt "  Base measure     : N(" as res %8.4g e(mu0) as txt ", " ///
        as res %8.4g e(s0sq) as txt ")"
    di as txt "  Sampler          : " as res e(chains) as txt " chain(s), " ///
        as res e(iters) as txt " sweeps, " as res e(burnin) as txt " burn-in"
    di as txt "  Assignment moves : " as res "`marginal'" as txt " covariance"
    di as txt "  Retained draws   : " as res e(ndraws)
    di as txt ""
    di as txt "  Posterior E[# groups] : " as res %6.3f e(m_mean)
    di as txt "{hline 72}"
    matlist e(cells), format(%9.4f) twidth(14)

    if ("`coclust'" == "") {
        di ""
        di as txt "  Posterior co-clustering probabilities:"
        matlist e(coclust), format(%6.2f) twidth(12)
    }
    di ""
    di as txt "{p 2 2 2}Intervals are posterior percentiles and already " ///
        "marginalize over the partition. A diffuse co-clustering matrix is " ///
        "not a failure: with few, correlated cells it is the correct report " ///
        "that the fine grouping is genuinely uncertain.{p_end}"
end
