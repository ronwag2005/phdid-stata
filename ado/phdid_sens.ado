*! version 0.1.0  phdid sensitivity -- regularisation paths
*! Rohan Wagle and Parush Arora
program define phdid_sens, rclass
    version 17.0
    syntax [, ALPHA(numlist >0 sort) LAMBDA(numlist >=0 sort) BY(string) ///
              TYPE(string) ITERs(integer 4000) BURNin(integer 1000) ///
              THIN(integer 1) CHAINs(integer 1) MARGinal(string) ///
              MU0(string) S0SQ(string) INIT(string) Level(cilevel) ///
              SEED(string) noGRaph GRAPHopts(string asis) SAVing(string) * ]

    // twoway options may be given directly as well as through
    // graphopts(), so that name() and title() work where you expect
    local graphopts `"`graphopts' `options'"'

    if ("`type'" == "") local type overall
    if !inlist("`type'", "overall", "cells") {
        di as err "type({bf:`type'}) not recognised; choose {bf:overall} or {bf:cells}"
        exit 198
    }
    if ("`by'" == "") local by lambda
    if !inlist("`by'", "lambda", "m") {
        di as err "by({bf:`by'}) not recognised; choose {bf:lambda} or {bf:m}"
        exit 198
    }
    if ("`marginal'" == "") local marginal exact
    if ("`init'" == "") local init singletons
    if ("`seed'" != "") set seed `seed'

    // Which path?  An explicit alpha() asks for the Dirichlet process
    // path; an explicit lambda() or by(m) asks for the l0 path; with
    // neither, follow whichever estimator was last run.
    local which ""
    if ("`alpha'" != "") local which alpha
    if ("`lambda'" != "" | "`by'" == "m") local which lambda
    if ("`alpha'" != "" & ("`lambda'" != "")) {
        di as err "specify {bf:alpha()} or {bf:lambda()}, not both"
        exit 198
    }
    if ("`which'" == "") {
        if ("`e(cmd)'" == "phdid_bayes") local which alpha
        else local which lambda
    }

    if ("`which'" == "alpha") {
        if ("`alpha'" == "") local alpha 0.25 0.5 1 2 4 8
        local mu0v = cond("`mu0'" == "", ".", "`mu0'")
        local s0v  = cond("`s0sq'" == "", ".", "`s0sq'")
        local exact = ("`marginal'" == "exact")
        local np : word count `alpha'
        di as txt "running `np' chains across alpha ..."
        mata: phdid_alphapath_run((`=subinstr("`alpha'", " ", ",", .)'), ///
            `mu0v', `s0v', `iters', `burnin', `thin', `chains', `exact', ///
            "`init'", `level', "`type'")
        local par alpha
        local parlab "concentration parameter alpha"
        local fixed "DP prior; `iters' iterations, `burnin' burn-in, `marginal' marginal"
    }
    else {
        if ("`lambda'" == "" & "`by'" == "lambda") local lambda 0 0.05 0.1 0.25 0.5 1 2 4
        if ("`by'" == "m") local lambda 0
        mata: phdid_lambdapath_run((`=subinstr("`lambda'", " ", ",", .)'), ///
            `level', "`type'", "`by'")
        local par lambda
        local parlab "l0 penalty lambda"
        if ("`by'" == "m") local parlab "number of groups m"
        local fixed "l0 estimator, step-2 GLS at each point"
    }

    tempname path
    matrix `path' = r(path)
    local np = r(npoints)

    di ""
    di as txt "{hline 72}"
    di as txt "Regularisation path over the `parlab'"
    di as txt "{hline 72}"
    di as txt "  `fixed'"
    di as txt ""
    matlist `path', format(%9.4f) twidth(10)
    di as txt "{hline 72}"
    if ("`type'" == "overall") {
        di as txt "{p 2 2 2}The aggregate is the quantity least sensitive " ///
            "to over-pooling, so a flat path here is weak evidence.  Run " ///
            "{bf:type(cells)} for the per-cell view, which is where " ///
            "over-pooling shows up.{p_end}"
    }
    else {
        di as txt "{p 2 2 2}One column per cohort-time effect.  Cells that " ///
            "move a long way across the path are the ones whose reported " ///
            "value is doing the work of the prior rather than of the " ///
            "data.{p_end}"
    }
    if ("`which'" == "lambda" & "`by'" == "lambda") {
        di as txt "{p 2 2 2}Note: the Appendix B stopping rule compares " ///
            "each merge cost against lambda times |A||B|, and that factor " ///
            "grows as groups absorb one another, so the effective " ///
            "threshold is not monotone in the merge order and some group " ///
            "counts are unreachable by any single lambda.  {bf:by(m)} " ///
            "walks the path one merge at a time instead.{p_end}"
    }

    if ("`saving'" != "") {
        preserve
        phdid_path2data `path'
        qui save `saving'
        restore
        di as txt "path saved to " as res "`saving'"
    }
    if ("`graph'" == "") {
        phdid_graph_path `path', param(`par') type(`type') ///
            parlab("`parlab'") level(`level') by(`by') `graphopts'
    }

    return matrix path = `path'
    return local parameter "`par'"
    return local path_type "`type'"
    return scalar npoints = `np'
end

// Turn a returned path matrix into the dataset in memory, so that it can
// be saved or graphed with ordinary twoway commands.  Column names come
// from the matrix, qualified by the equation name when there is one and
// uniquified when two cells still collide -- the cohort-time labels are
// not themselves legal variable names.
program define phdid_path2data
    version 17.0
    args M
    local nr = rowsof(`M')
    local nc = colsof(`M')
    local cn : colnames `M'
    local ce : coleq `M'
    clear
    qui set obs `nr'
    local vars ""
    local labs ""
    forvalues j = 1/`nc' {
        local nm : word `j' of `cn'
        local eq : word `j' of `ce'
        local full "`nm'"
        if ("`eq'" != "" & "`eq'" != "_") local full "`eq':`nm'"
        local v = strtoname("`full'")
        local taken : list posof "`v'" in vars
        if (`taken' > 0) local v "`v'_`j'"
        qui gen double `v' = .
        forvalues i = 1/`nr' {
            qui replace `v' = `M'[`i', `j'] in `i'
        }
        label variable `v' "`full'"
        local vars `vars' `v'
        local labs `"`labs' `"`full'"'"'
    }
    c_local phdid_pvars `vars'
    c_local phdid_plabs `"`labs'"'
end

program define phdid_graph_path
    version 17.0
    syntax anything(name=M) , PARam(string) TYPE(string) PARLab(string) ///
        Level(cilevel) BY(string) [ * ]
    local opts `"`options'"'

    preserve
    phdid_path2data `M'
    local pv `phdid_pvars'
    local pl `"`phdid_plabs'"'
    local xv : word 1 of `pv'
    local mv : word 2 of `pv'
    local xtra ""
    // by(m) walks the agglomeration path one merge at a time, so the
    // group count is the axis and the implied lambda is incidental.
    if ("`by'" == "m") {
        local xv `mv'
        qui su `mv'
        local xtra xlabel(`r(min)'(1)`r(max)')
    }

    if ("`type'" == "overall") {
        local av : word 3 of `pv'
        local lv : word 5 of `pv'
        local uv : word 6 of `pv'
        twoway (rarea `uv' `lv' `xv', color(gs13)) ///
               (line `av' `xv', lcolor(navy) lwidth(medthick)) ///
               (line `mv' `xv', lcolor(maroon) lpattern(dash) yaxis(2)), ///
            legend(order(2 "overall effect" 1 "`level'% interval" ///
                3 "groups (right axis)") rows(1) size(small)) ///
            xtitle("`parlab'") ytitle("overall effect") ///
            ytitle("number of groups", axis(2)) ///
            yline(0, lcolor(gs8) lpattern(dot)) ///
            title("phdid: sensitivity of the aggregate") ///
            note("A flat line here is weak evidence; see type(cells).") ///
            `xtra' graphregion(color(white)) `opts'
    }
    else {
        local pal navy maroon forest_green dkorange teal cranberry ///
            lavender olive_teal sienna purple gold eltblue magenta gs6
        local nv : word count `pv'
        local lines ""
        local lg ""
        forvalues j = 3/`nv' {
            local v : word `j' of `pv'
            local L : word `j' of `pl'
            local kk = `j' - 2
            local ci = mod(`kk' - 1, 14) + 1
            local col : word `ci' of `pal'
            local lines `lines' (line `v' `xv', lcolor(`col') lwidth(medium))
            local lg `"`lg' `kk' "`L'""'
        }
        twoway `lines', legend(order(`lg') size(vsmall) cols(4)) ///
            xtitle("`parlab'") ytitle("cell effect") ///
            yline(0, lcolor(gs8) lpattern(dot)) ///
            title("phdid: sensitivity of each cohort-time effect") ///
            `xtra' graphregion(color(white)) `opts'
    }
    restore
end
