*! version 0.1.0  phdid fit -- restricted GLS under a given partition
*! Rohan Wagle and Parush Arora
program define phdid_fit, eclass
    version 17.0
    syntax , PARTition(numlist integer >0) [ Level(cilevel) noTABle ]

    local K : word count `partition'
    tempname zm
    matrix `zm' = J(`K', 1, .)
    local i 0
    foreach z of local partition {
        local ++i
        matrix `zm'[`i', 1] = `z'
    }

    mata: phdid_fit_run(st_matrix("`zm'"))

    tempname bmat Vmat cmat phimat pz
    matrix `bmat'   = r(b)
    matrix `Vmat'   = r(V)
    matrix `cmat'   = r(cells)
    matrix `phimat' = r(phi)
    matrix `pz'     = r(partition)
    local m    = r(m)
    local dev  = r(deviance)
    local df   = r(df)
    local pv   = r(p)
    local att  = r(att)
    local ase  = r(att_se)
    local Kr   = r(K)
    local pstr "`r(partition_str)'"

    ereturn post `bmat' `Vmat', depname("CATT")
    ereturn local cmd      "phdid_fit"
    ereturn local cmdline  "phdid fit, partition(`partition')"
    ereturn local title    "restricted GLS under a supplied partition"
    ereturn local partition_str "`pstr'"
    ereturn matrix cells = `cmat'
    ereturn matrix phi   = `phimat'
    ereturn matrix partition = `pz'
    ereturn scalar m        = `m'
    ereturn scalar deviance = `dev'
    ereturn scalar df       = `df'
    ereturn scalar p        = `pv'
    ereturn scalar att      = `att'
    ereturn scalar att_se   = `ase'
    ereturn scalar K        = `Kr'

    di ""
    di as txt "{hline 72}"
    di as txt "Restricted GLS under a supplied partition"
    di as txt "{hline 72}"
    di as txt "  Cells (K)   : " as res `Kr'
    di as txt "  Groups (m)  : " as res `m'
    di as txt "  Partition   : " as res "`pstr'"
    di as txt "{hline 72}"
    if ("`table'" == "") {
        matlist e(cells), format(%9.4f) twidth(14)
    }
    di as txt ""
    di as txt "  Overall ATT : " as res %9.4f `att' ///
        as txt "  (se " as res %7.4f `ase' as txt ")"
    di as txt "  Deviance    : " as res %9.4f `dev' ///
        as txt "  chi2(" as res `df' as txt ") against the flexible fit," ///
        as txt " p = " as res %6.4f `pv'
    di as txt "{hline 72}"
    di as txt "{p 2 2 2}The deviance is the cost of imposing this " ///
        "partition relative to estimating every cell freely.  A small " ///
        "p-value says the restriction is rejected by the data.  Because " ///
        "the partition was not selected here, that test is an honest one; " ///
        "the same statistic after {bf:phdid l0} is not, since the " ///
        "partition was chosen to make it small.{p_end}"
end
