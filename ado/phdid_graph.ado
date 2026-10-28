*! version 0.1.0  phdid graph -- co-clustering, grouped effects, BIC path
*! Rohan Wagle and Parush Arora
program define phdid_graph
    version 17.0
    gettoken kind 0 : 0, parse(" ,")
    if ("`kind'" == "" | substr("`kind'", 1, 1) == ",") {
        local 0 `"`kind'`0'"'
        if ("`e(cmd)'" == "phdid_bayes") local kind coclust
        else local kind groups
    }
    if !inlist("`kind'", "coclust", "groups", "bic") {
        di as err "{bf:phdid graph `kind'} not recognised"
        di as txt "{p 4 4 2}Choose {bf:coclust} (posterior co-clustering " ///
            "matrix), {bf:groups} (flexible against grouped effects) or " ///
            "{bf:bic} (the l0 selection path).{p_end}"
        exit 198
    }
    if ("`kind'" == "coclust") phdid_graph_cc `0'
    if ("`kind'" == "groups")  phdid_graph_gr `0'
    if ("`kind'" == "bic")     phdid_graph_bic `0'
end

// The co-clustering matrix, drawn as a grid of shaded squares.  Stata's
// contour plots interpolate between grid points, which would imply a
// continuum that does not exist here, so each cell is drawn exactly as
// one marker and the shading is banded.
program define phdid_graph_cc
    version 17.0
    syntax [, noVALues * ]

    tempname CC
    if ("`e(cmd)'" == "phdid_bayes") {
        matrix `CC' = e(coclust)
    }
    else {
        capture matrix `CC' = r(coclust)
    }
    capture local K = rowsof(`CC')
    if (_rc | "`K'" == "" | `K' == 0) {
        di as err "no co-clustering matrix available"
        di as txt "{p 4 4 2}Run {bf:phdid bayes} (or {bf:phdid enumerate}) " ///
            "first.{p_end}"
        exit 301
    }
    local rn : rownames `CC'
    local re : roweq `CC'

    preserve
    clear
    qui set obs `=`K' * `K''
    qui gen double _row = .
    qui gen double _col = .
    qui gen double _p   = .
    qui gen str4   _lab = ""
    local n 0
    forvalues a = 1/`K' {
        forvalues b = 1/`K' {
            local ++n
            qui replace _row = `a' in `n'
            qui replace _col = `b' in `n'
            qui replace _p   = `CC'[`a', `b'] in `n'
        }
    }
    qui replace _lab = string(_p, "%3.2f")

    local xlab ""
    forvalues a = 1/`K' {
        local nm : word `a' of `rn'
        local eq : word `a' of `re'
        local full "`nm'"
        if ("`eq'" != "" & "`eq'" != "_") local full "`eq':`nm'"
        local xlab `xlab' `a' "`full'"
    }

    // Each cell is one bar, one data unit wide and one data unit tall,
    // so the grid is drawn exactly.  Contour plots interpolate between
    // grid points, which would imply a continuum between cohort-time
    // cells that does not exist.
    qui gen double _ylo = _row - 0.5
    qui gen double _yhi = _row + 0.5
    local pal gs15 gs12 gs9 gs6 gs3
    local lay ""
    local lg ""
    forvalues q = 1/5 {
        local lo = (`q' - 1) / 5
        local hi = `q' / 5
        local cnd "_p > `lo' & _p <= `hi'"
        if (`q' == 1) local cnd "_p >= 0 & _p <= `hi'"
        local col : word `q' of `pal'
        local lay `lay' (rbar _ylo _yhi _col if `cnd', barwidth(1) ///
            color(`col') lcolor(white) lwidth(vthin))
        local lg `lg' `q' "`:di %3.1f `lo''-`:di %3.1f `hi''"
    }
    if ("`values'" == "") {
        local lay `lay' (scatter _row _col if _p <= 0.6, msymbol(none) ///
            mlabel(_lab) mlabsize(small) mlabposition(0) mlabcolor(black))
        local lay `lay' (scatter _row _col if _p > 0.6, msymbol(none) ///
            mlabel(_lab) mlabsize(small) mlabposition(0) mlabcolor(white))
    }

    twoway `lay', ///
        legend(order(`lg') title("Pr(same group)", size(small)) ///
            rows(1) size(vsmall) position(6) region(lstyle(none))) ///
        xlabel(`xlab', angle(45) labsize(vsmall) notick) ///
        ylabel(`xlab', angle(0) labsize(vsmall) notick) yscale(reverse) ///
        xscale(range(0.5 `=`K' + 0.5') lstyle(none)) ///
        yscale(range(0.5 `=`K' + 0.5') lstyle(none)) ///
        xtitle("") ytitle("") aspectratio(1) ///
        title("phdid: posterior co-clustering probabilities") ///
        note("Darker cells are ones the posterior keeps together.") ///
        graphregion(color(white)) `options'
    restore
end

// Flexible against grouped estimates, cell by cell.  The point of the
// estimator is visible here: within a recovered group the grouped
// estimates collapse onto one line with a shorter interval.
program define phdid_graph_gr
    version 17.0
    syntax [, Level(cilevel) * ]

    tempname CM
    if !inlist("`e(cmd)'", "phdid_l0", "phdid_bayes", "phdid_fit", ///
               "phdid_partition") {
        di as err "no phdid estimation results found"
        exit 301
    }
    matrix `CM' = e(cells)
    local K = rowsof(`CM')
    local rn : rownames `CM'
    local cn : colnames `CM'
    local isl0 = ("`e(cmd)'" != "phdid_bayes")
    local zq = invnormal(1 - (1 - `level' / 100) / 2)

    preserve
    clear
    qui set obs `K'
    qui gen double _k   = _n
    qui gen double _kf  = _n - 0.14
    qui gen double _kg  = _n + 0.14
    qui gen double _fl  = .
    qui gen double _fse = .
    qui gen double _gr  = .
    qui gen double _gse = .
    qui gen double _glo = .
    qui gen double _ghi = .
    qui gen int    _grp = .
    forvalues i = 1/`K' {
        qui replace _fl  = `CM'[`i', 1] in `i'
        if (`isl0') {
            qui replace _fse = `CM'[`i', 2] in `i'
            qui replace _grp = `CM'[`i', 3] in `i'
            qui replace _gr  = `CM'[`i', 4] in `i'
            qui replace _gse = `CM'[`i', 5] in `i'
        }
        else {
            qui replace _gr  = `CM'[`i', 2] in `i'
            qui replace _gse = `CM'[`i', 3] in `i'
            qui replace _glo = `CM'[`i', 4] in `i'
            qui replace _ghi = `CM'[`i', 5] in `i'
        }
    }
    if (`isl0') {
        qui replace _glo = _gr - `zq' * _gse
        qui replace _ghi = _gr + `zq' * _gse
    }
    qui gen double _flo = _fl - `zq' * _fse
    qui gen double _fhi = _fl + `zq' * _fse

    local re : roweq `CM'
    local xlab ""
    forvalues i = 1/`K' {
        local nm : word `i' of `rn'
        local eq : word `i' of `re'
        local full "`nm'"
        if ("`eq'" != "" & "`eq'" != "_") local full "`eq':`nm'"
        local xlab `xlab' `i' "`full'"
    }
    local lab2 = cond(`isl0', "grouped", "posterior mean")

    local flex ""
    if (`isl0') {
        local flex (rcap _fhi _flo _kf, lcolor(gs9) msize(small)) ///
                   (scatter _fl _kf, msymbol(oh) mcolor(gs6))
        local lgf 1 "flexible, `level'% CI"
        local n1 2
    }
    else {
        local flex (scatter _fl _kf, msymbol(oh) mcolor(gs6))
        local lgf ""
        local n1 1
    }

    twoway `flex' ///
           (rcap _ghi _glo _kg, lcolor(navy) msize(small)) ///
           (scatter _gr _kg, msymbol(O) mcolor(navy)), ///
        legend(order(`=`n1'' "flexible" `=`n1'+2' "`lab2'") ///
            rows(1) size(small)) ///
        xlabel(`xlab', angle(45) labsize(small)) xtitle("") ///
        ytitle("effect on the outcome") ///
        yline(0, lcolor(gs8) lpattern(dot)) ///
        xscale(range(0.4 `=`K' + 0.6')) ///
        title("phdid: flexible against grouped cohort-time effects") ///
        note("Grey: the unrestricted first-stage effects. " ///
             "Navy: after partial pooling.") ///
        graphregion(color(white)) `options'
    restore
end

// The l0 selection path: deviance and BIC against the number of groups.
program define phdid_graph_bic
    version 17.0
    syntax [, * ]

    if ("`e(cmd)'" != "phdid_l0") {
        di as err "no phdid l0 results found"
        exit 301
    }
    tempname PM
    matrix `PM' = e(path)
    local nr = rowsof(`PM')
    local msel = e(m)

    preserve
    clear
    qui set obs `nr'
    qui gen int    _m   = .
    qui gen double _dev = .
    qui gen double _bic = .
    forvalues i = 1/`nr' {
        qui replace _m   = `PM'[`i', 1] in `i'
        qui replace _dev = `PM'[`i', 2] in `i'
        qui replace _bic = `PM'[`i', 3] in `i'
    }

    twoway (connected _bic _m, lcolor(navy) mcolor(navy) msymbol(O)) ///
           (connected _dev _m, lcolor(maroon) mcolor(maroon) ///
                msymbol(Oh) lpattern(dash)), ///
        xline(`msel', lcolor(gs8) lpattern(dot)) ///
        legend(order(1 "BIC" 2 "deviance") rows(1) size(small)) ///
        xtitle("number of groups m") ytitle("") ///
        xlabel(1(1)`nr') ///
        title("phdid: the l0 selection path") ///
        note("Vertical line: m = `msel', the selected partition. " ///
             "BIC penalty uses n = `e(nbic)'.") ///
        graphregion(color(white)) `options'
    restore
end
