*! version 0.1.0  phdid l0 -- l0-penalised partial-homogeneity estimator
*! Rohan Wagle and Parush Arora
program define phdid_l0, eclass
    version 17.0
    syntax [, SELect(string) Lambda(real -1) M(integer -1) NBIC(real -1) ///
              Level(cilevel) noPATH ]

    if ("`select'" == "") local select "bic"
    if (!inlist("`select'", "bic", "lambda", "m")) {
        di as err "select() must be {bf:bic}, {bf:lambda} or {bf:m}"
        exit 198
    }
    if ("`select'" == "lambda" & `lambda' < 0) {
        di as err "select(lambda) requires a non-negative {bf:lambda()}"
        exit 198
    }
    if ("`select'" == "m" & `m' < 1) {
        di as err "select(m) requires {bf:m()} of at least 1"
        exit 198
    }

    mata: st_numscalar("r(K)", phdid_K())
    local K = r(K)
    if (`nbic' < 0) local nbic = `K'
    if ("`select'" == "m" & `m' > `K') {
        di as err "m() may not exceed K = `K'"
        exit 198
    }

    mata: phdid_l0_run("`select'", `lambda', `m', `nbic')

    tempname bmat Vmat cmat pmat phimat zmat
    matrix `bmat'   = r(b)
    matrix `Vmat'   = r(V)
    matrix `cmat'   = r(cells)
    matrix `pmat'   = r(path)
    matrix `phimat' = r(phi)
    matrix `zmat'   = r(partition)
    local msel     = r(m)
    local dev      = r(deviance)
    local bicv     = r(bic)
    local att      = r(att)
    local attse    = r(att_se)
    local pstr     "`r(partition_str)'"

    ereturn post `bmat' `Vmat', depname("CATT")
    ereturn local cmd      "phdid_l0"
    ereturn local cmdline  "phdid l0 `0'"
    ereturn local title    "l0-penalised partial-homogeneity estimator"
    ereturn local select   "`select'"
    ereturn local partition_str "`pstr'"
    ereturn matrix cells = `cmat'
    ereturn matrix path  = `pmat'
    ereturn matrix phi   = `phimat'
    ereturn matrix partition = `zmat'
    ereturn scalar m        = `msel'
    ereturn scalar deviance = `dev'
    ereturn scalar bic      = `bicv'
    ereturn scalar att      = `att'
    ereturn scalar att_se   = `attse'
    ereturn scalar K        = `K'
    ereturn scalar nbic     = `nbic'
    if ("`select'" == "lambda") ereturn scalar lambda = `lambda'

    phdid_l0_display, level(`level') `nopath'
end

program define phdid_l0_display
    version 17.0
    syntax [, Level(cilevel) noPATH ]

    di ""
    di as txt "{hline 72}"
    di as txt "l0-penalised partial-homogeneity estimator"
    di as txt "{hline 72}"
    di as txt "  Cells (K)   : " as res e(K)
    if ("`e(select)'" == "bic") {
        di as txt "  Groups (m)  : " as res e(m) as txt "   (selected by BIC, n = " as res e(nbic) as txt " in the penalty)"
    }
    else if ("`e(select)'" == "lambda") {
        di as txt "  Groups (m)  : " as res e(m) as txt "   (Appendix B rule at lambda = " as res e(lambda) as txt ")"
    }
    else {
        di as txt "  Groups (m)  : " as res e(m) as txt "   (fixed by hand)"
    }
    di as txt "  Partition   : " as res "`e(partition_str)'"
    di as txt "{hline 72}"
    matlist e(cells), format(%9.4f) twidth(14)
    di ""
    di as txt "  var_ratio is the grouped sampling variance relative to flexible;"
    di as txt "  below 1 means pooling bought precision for that cell."
    di ""
    di as txt "  Overall ATT : " as res %9.4f e(att) as txt "  (se " as res %7.4f e(att_se) as txt ")"

    if ("`nopath'" == "") {
        di ""
        di as txt "  Agglomeration path:"
        matlist e(path), format(%12.5g) twidth(6)
    }
    di ""
    di as txt "{p 2 2 2}These intervals condition on the selected partition " ///
        "being correct. They can under-cover, badly so when the partition is " ///
        "uncertain. Use {bf:phdid bayes} for inference that averages over " ///
        "partitions.{p_end}"
end
