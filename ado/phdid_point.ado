*! version 0.1.0  phdid partition -- point partition from the posterior
*! Rohan Wagle and Parush Arora
program define phdid_point, rclass
    version 17.0
    syntax [, LOSS(string) POST ]

    if ("`loss'" == "") local loss vi
    if !inlist("`loss'", "vi", "binder") {
        di as err "loss({bf:`loss'}) not recognised; choose {bf:vi} or {bf:binder}"
        exit 198
    }
    if ("`e(cmd)'" != "phdid_bayes") {
        di as err "no phdid bayes results found"
        exit 301
    }

    mata: phdid_point_run("`loss'")

    tempname zz bb VV ccm
    matrix `zz'  = r(partition)
    matrix `bb'  = r(b)
    matrix `VV'  = r(V)
    matrix `ccm' = r(coclust)
    local m    = r(m)
    local lv   = r(loss_value)
    local nu   = r(n_unique)
    local pp   = r(post_prob)
    local att  = r(att)
    local pstr "`r(partition_str)'"

    local lname = cond("`loss'" == "vi", "variation of information", "Binder")

    di ""
    di as txt "{hline 72}"
    di as txt "Point partition minimising expected `lname' loss"
    di as txt "{hline 72}"
    di as txt "  Groups (m)                 : " as res `m'
    di as txt "  Distinct partitions sampled: " as res %9.0fc `nu'
    di as txt "  Posterior mass on this one : " as res %6.4f `pp'
    di as txt "  Expected loss              : " as res %9.5f `lv'
    di as txt "  Overall effect under it    : " as res %9.5f `att'
    di as txt ""
    di as txt "  `pstr'"
    di as txt "{hline 72}"
    di as txt "{p 2 2 2}A single partition is a summary of the posterior, " ///
        "not a finding: the mass on any one grouping is small, which is " ///
        "why the interval from {bf:phdid aggregate} should still be the " ///
        "one reported.  The co-clustering matrix says which parts of this " ///
        "grouping are firm.{p_end}"

    return matrix partition = `zz'
    return matrix coclust = `ccm'
    return scalar m = `m'
    return scalar loss_value = `lv'
    return scalar post_prob = `pp'
    return local loss "`loss'"
    return local partition_str "`pstr'"

    if ("`post'" != "") {
        ereturn post `bb' `VV', depname("CATT")
        ereturn local cmd     "phdid_partition"
        ereturn local cmdline "phdid partition, loss(`loss')"
        ereturn local partition_str "`pstr'"
        ereturn scalar m = `m'
        di ""
        di as txt "{p 2 2 2}Restricted GLS fit under this partition posted " ///
            "to {bf:e(b)} and {bf:e(V)}.  These are conditional standard " ///
            "errors.{p_end}"
    }
    else {
        return matrix b = `bb'
        return matrix V = `VV'
    }
end
