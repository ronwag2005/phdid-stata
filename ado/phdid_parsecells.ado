*! version 0.1.0  phdid_parsecells -- recover (g,t) from coefficient names
*! Rohan Wagle and Parush Arora
*
* First-stage commands name their cohort-time coefficients differently.
* This recognises the common conventions and returns the indices of the
* post-treatment cells (t >= g) with the recovered g and t, and where the
* first stage also posts cohort weights, those too.
*
* Recognised:
*   csdid    eq "g2004", name "t_2003_2004"   -> g from eq, t from last number
*            plus an eq "wgt" carrying cohort shares as "w2004_2004"
*   generic  name "2004:2004" or "2004_2004", or a numeric equation
*            name carrying the cohort and the column name the period
*   factor   name "2004.g#2005.t"
*
* When nothing matches, every column is kept and (g,t) are left empty,
* which disables only the aggregations that need them.
program define phdid_parsecells, rclass
    version 17.0
    syntax , LABels(string) [ EQnames(string) ]

    local keep  ""
    local glist ""
    local tlist ""
    local nmout ""
    local wix   ""
    local wgt   ""
    local matched 0
    local n : word count `labels'

    * ---------- pass 1: csdid-style, g from the equation name ----------
    if (`"`eqnames'"' != "") {
        local anyg 0
        local i 0
        foreach nm of local labels {
            local ++i
            local eq : word `i' of `eqnames'
            if regexm("`eq'", "^g([0-9]+)$") {
                local anyg 1
                local g = real(regexs(1))
                * t is the final run of digits in the coefficient name
                local t = .
                if regexm("`nm'", "([0-9]+)$") {
                    local t = real(regexs(1))
                }
                if (`t' < . & `t' >= `g') {
                    local keep  "`keep' `i'"
                    local glist "`glist' `g'"
                    local tlist "`tlist' `t'"
                    local nmout "`nmout' `g':`t'"
                }
            }
            else if regexm("`eq'", "^([0-9]+)$") {
                * a bare numeric equation name is the cohort: this is how
                * a hand-built matrix with colnames "2004:2005" arrives,
                * Stata reading the colon as the equation separator
                local g = real(regexs(1))
                local t = .
                if regexm("`nm'", "([0-9]+)$") {
                    local t = real(regexs(1))
                }
                if (`t' < . & `t' >= `g') {
                    local anyg 1
                    local keep  "`keep' `i'"
                    local glist "`glist' `g'"
                    local tlist "`tlist' `t'"
                    local nmout "`nmout' `g':`t'"
                }
            }
            else if ("`eq'" == "wgt") {
                if regexm("`nm'", "^w([0-9]+)_([0-9]+)$") {
                    local wix "`wix' `i'"
                    local wgt "`wgt' `=regexs(1)':`=regexs(2)'"
                }
            }
        }
        if (`anyg') local matched 1
    }

    * ---------- pass 2: name carries both g and t ----------
    if (`matched' == 0) {
        local keep "" 
        local glist ""
        local tlist ""
        local nmout ""
        local ok 1
        local i 0
        foreach nm of local labels {
            local ++i
            local g .
            local t .
            if regexm("`nm'", "^([0-9]+)[:_]([0-9]+)$") {
                local g = real(regexs(1))
                local t = real(regexs(2))
            }
            else if regexm("`nm'", "([0-9]+)\.[a-zA-Z_]+#([0-9]+)\.[a-zA-Z_]+") {
                local g = real(regexs(1))
                local t = real(regexs(2))
            }
            else {
                local ok 0
                continue, break
            }
            if (`t' >= `g') {
                local keep  "`keep' `i'"
                local glist "`glist' `g'"
                local tlist "`tlist' `t'"
                local nmout "`nmout' `g':`t'"
            }
        }
        if (`ok') local matched 1
    }

    * ---------- fallback: keep everything, no (g,t) ----------
    if (`matched' == 0) {
        local keep ""
        local nmout ""
        forvalues j = 1/`n' {
            local keep "`keep' `j'"
            local nm : word `j' of `labels'
            local nmout "`nmout' `nm'"
        }
        local glist ""
        local tlist ""
    }

    return local keep   "`keep'"
    return local g      "`glist'"
    return local t      "`tlist'"
    return local names  "`nmout'"
    return local wix    "`wix'"
    return local wgtkey "`wgt'"
    return scalar matched = `matched'
end
