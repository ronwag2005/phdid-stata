*! version 0.1.0  phdid setup -- assemble the first-stage input (Lemma 1)
*! Rohan Wagle and Parush Arora
program define phdid_setup, rclass
    version 17.0
    syntax [anything] [, B(name) V(name) CELLS(string) Weights(name) ///
        NUnits(real -1) NOWARN ]

    tempname bb vv ww gg tt

    * ---- source the first stage ----
    if ("`b'" != "" | "`v'" != "") {
        if ("`b'" == "" | "`v'" == "") {
            di as err "specify both {bf:b()} and {bf:v()}, or neither"
            exit 198
        }
        matrix `bb' = `b'
        matrix `vv' = `v'
        local src "supplied matrices"
    }
    else {
        if ("`e(cmd)'" == "") {
            di as err "no estimation results in memory"
            di as txt "{p 4 4 2}Run a heterogeneity-robust first stage such " ///
                "as {bf:csdid} or {bf:jwdid} first, or supply the estimates " ///
                "with {bf:b()} and {bf:v()}.{p_end}"
            exit 301
        }
        matrix `bb' = e(b)
        matrix `vv' = e(V)
        local src "e(b) and e(V) from {bf:`e(cmd)'}"
        if (`nunits' == -1) {
            cap confirm scalar e(N)
            if (_rc == 0) local nunits = e(N)
        }
    }

    * ---- cell labels, and (g,t) where they can be recovered ----
    local K = colsof(`bb')
    local labs : colnames `bb'
    local eqs  : coleq `bb'

    * drop pre-treatment cells, and pick up cohort weights where the
    * first stage posts them
    phdid_parsecells, labels(`labs') eqnames(`eqs')
    local keep   "`r(keep)'"
    local glist  "`r(g)'"
    local tlist  "`r(t)'"
    local pnames "`r(names)'"
    local wix    "`r(wix)'"
    local wgtkey "`r(wgtkey)'"
    local nkeep : word count `keep'

    * cells() overrides whatever the coefficient names said.  The names
    * supplied are parsed in turn, so "2004:2005" still yields a cohort
    * and a period and still drops the pre-treatment cells.
    if ("`cells'" != "") {
        local labs "`cells'"
        local nlab : word count `labs'
        if (`nlab' != `K') {
            di as err "cells() has `nlab' names but b() has `K' columns"
            exit 198
        }
        phdid_parsecells, labels(`labs')
        local keep   "`r(keep)'"
        local glist  "`r(g)'"
        local tlist  "`r(t)'"
        local pnames "`r(names)'"
        local nkeep : word count `keep'
        local wix    ""
        local wgtkey ""
    }

    if (`nkeep' < 2) {
        di as err "fewer than two cohort-time cells were identified"
        exit 198
    }

    * ---- subset ----
    tempname bs vs
    matrix `bs' = J(1, `nkeep', 0)
    matrix `vs' = J(`nkeep', `nkeep', 0)
    local newlab ""
    local j 0
    foreach i of local keep {
        local ++j
        matrix `bs'[1, `j'] = `bb'[1, `i']
        local nm : word `j' of `pnames'
        local newlab "`newlab' `nm'"
        local l 0
        foreach ii of local keep {
            local ++l
            matrix `vs'[`j', `l'] = `vv'[`i', `ii']
        }
    }
    matrix colnames `bs' = `newlab'

    * ---- weights ----
    * If the first stage posted cohort shares (csdid does), use them: the
    * overall ATT is then the population-weighted average, matching the
    * paper's first application. Otherwise weight the cells equally.
    if ("`weights'" != "") {
        matrix `ww' = `weights'
        if (rowsof(`ww') > 1 & colsof(`ww') == 1) matrix `ww' = `ww''
        if (rowsof(`ww') != 1 | colsof(`ww') != `nkeep') {
            di as err "weights() must have one element per cell"
            di as txt "{p 4 4 2}There are `nkeep' post-treatment cells, " ///
                "but {bf:weights()} is " as res "`=rowsof(`ww')'" ///
                as txt " by " as res "`=colsof(`ww')'" as txt ".{p_end}"
            exit 503
        }
        local wsrc "supplied"
    }
    else if ("`wix'" != "" & "`wgtkey'" != "") {
        matrix `ww' = J(1, `nkeep', .)
        local j 0
        foreach nm of local pnames {
            local ++j
            local pos 0
            local l 0
            foreach kk of local wgtkey {
                local ++l
                if ("`kk'" == "`nm'") local pos = `l'
            }
            if (`pos' > 0) {
                local col : word `pos' of `wix'
                matrix `ww'[1, `j'] = `bb'[1, `col']
            }
        }
        cap assert !missing(`ww'[1,1])
        if (_rc) {
            matrix `ww' = J(1, `nkeep', 1)
            local wsrc "equal"
        }
        else {
            local wsrc "cohort shares from the first stage"
        }
    }
    else {
        matrix `ww' = J(1, `nkeep', 1)
        local wsrc "equal"
    }

    if ("`glist'" != "") {
        matrix `gg' = J(1, `nkeep', .)
        matrix `tt' = J(1, `nkeep', .)
        local j 0
        foreach gv of local glist {
            local ++j
            matrix `gg'[1, `j'] = `gv'
        }
        local j 0
        foreach tv of local tlist {
            local ++j
            matrix `tt'[1, `j'] = `tv'
        }
        local gname "`gg'"
        local tname "`tt'"
    }
    else {
        local gname ""
        local tname ""
    }

    global PHDID_LAB "`newlab'"
    mata: st_numscalar("r(K)", phdid_setup("`bs'", "`vs'", "PHDID_LAB", ///
        "`ww'", "`gname'", "`tname'", `nunits'))

    * ---- report ----
    mata: phdid_describe()
    local isdiag = r(isdiag)
    di ""
    di as txt "{hline 68}"
    di as txt "phdid: first-stage cohort-time effects"
    di as txt "{hline 68}"
    di as txt "  Source            : " as res "`src'"
    di as txt "  Cells (K)         : " as res r(K)
    di as txt "  Weights           : " as res "`wsrc'"
    if (`nunits' > 0) di as txt "  Units             : " as res `nunits'
    if (`isdiag') {
        di as txt "  Covariance        : " as res "diagonal"
    }
    else {
        di as txt "  Covariance        : " as res "full (cross-cell correlation retained)"
    }
    di as txt "  sd(tau)/median(se): " as res %6.2f r(snr_desc) as txt "   (descriptive)"
    di as txt "{hline 68}"
    matlist r(cells), format(%9.4f) twidth(14)

    if (`isdiag' & "`nowarn'" == "") {
        di ""
        di as txt "{p 2 2 2}Note: the covariance is diagonal, so the " ///
            "cohort-time effects are being treated as independent. Cells " ///
            "that share units, cohorts or calendar periods are generally " ///
            "correlated; ignoring this understates the posterior variance " ///
            "of aggregates such as the overall ATT and can misstate " ///
            "co-clustering probabilities. Supply the exact joint covariance " ///
            "where you can.{p_end}"
    }
    di ""
    di as txt "{p 2 2 2}Next: {bf:phdid test} asks whether there is " ///
        "heterogeneity to recover; {bf:phdid l0} and {bf:phdid bayes} " ///
        "recover it.{p_end}"

    return scalar K = r(K)
end
