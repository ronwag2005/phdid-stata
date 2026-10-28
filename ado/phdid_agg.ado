*! version 0.1.0  phdid aggregate -- overall, dynamic, group, calendar
*! Rohan Wagle and Parush Arora
program define phdid_agg, rclass
    version 17.0
    syntax [anything(name=type)] [, Level(cilevel) POST ]

    if ("`type'" == "") local type overall
    if !inlist("`type'", "overall", "dynamic", "group", "calendar") {
        di as err "aggregation {bf:`type'} not recognised"
        di as txt "{p 4 4 2}Choose {bf:overall}, {bf:dynamic} (by event " ///
            "time), {bf:group} (by cohort) or {bf:calendar} (by " ///
            "period).{p_end}"
        exit 198
    }

    if ("`e(cmd)'" == "phdid_bayes") {
        local src draws
        local note "partition averaged"
    }
    else if (inlist("`e(cmd)'", "phdid_l0", "phdid_fit", "phdid_partition")) {
        local src analytic
        local note "partition conditional"
    }
    else {
        di as err "no phdid estimation results found"
        di as txt "{p 4 4 2}Run {bf:phdid bayes}, {bf:phdid l0} or " ///
            "{bf:phdid fit} first.{p_end}"
        exit 301
    }

    mata: phdid_agg_run("`type'", `level', "`src'")

    tempname tbl bb VV AA
    matrix `tbl' = r(table)
    matrix `bb'  = r(b)
    matrix `VV'  = r(V)
    matrix `AA'  = r(A)
    local nterms = r(nterms)

    di ""
    di as txt "{hline 72}"
    if ("`type'" == "overall")  local tlab "Overall effect"
    if ("`type'" == "dynamic")  local tlab "Effects by event time"
    if ("`type'" == "group")    local tlab "Effects by cohort"
    if ("`type'" == "calendar") local tlab "Effects by period"
    di as txt "`tlab'" _col(44) as txt "`level'% CI, `note'"
    di as txt "{hline 72}"
    matlist `tbl', format(%9.5f) twidth(14)
    di as txt "{hline 72}"
    if ("`src'" == "draws") {
        di as txt "{p 2 2 2}Each term is a fixed linear map of the cell " ///
            "effects applied to every posterior draw, so the interval " ///
            "includes uncertainty about the grouping itself.{p_end}"
    }
    else {
        di as txt "{p 2 2 2}Standard errors condition on the selected " ///
            "partition being the true one.  When the partition is " ///
            "uncertain this understates the uncertainty; {bf:phdid bayes} " ///
            "marginalises over it instead.{p_end}"
    }

    return matrix table = `tbl'
    return matrix A = `AA'
    return local aggregation "`type'"
    return scalar nterms = `nterms'

    if ("`post'" != "") {
        local cl `e(cmdline)'
        ereturn post `bb' `VV', depname("ATT")
        ereturn local cmd      "phdid_aggregate"
        ereturn local cmdline  "phdid aggregate `type'"
        ereturn local parent   "`cl'"
        ereturn local aggregation "`type'"
        ereturn scalar level = `level'
        di ""
        di as txt "{p 2 2 2}Aggregates posted to {bf:e(b)} and {bf:e(V)}, " ///
            "replacing the cell-level fit; {bf:lincom} and {bf:test} now " ///
            "operate on the terms above.  Re-run the estimation command " ///
            "to get the cells back.{p_end}"
    }
    else {
        return matrix b = `bb'
        return matrix V = `VV'
    }
end
