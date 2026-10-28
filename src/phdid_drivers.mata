*! phdid_drivers.mata -- Stata-facing drivers for the phdid package
*! Rohan Wagle and Parush Arora
version 17.0

mata: mata set matastrict on

mata:

// =====================================================================
// INPUT ASSEMBLY  (paper Lemma 1)
//
// Every estimator in this package depends on the data only through
// R'Sigma^-1 tau and R'Sigma^-1 R.  Any consistent first-stage
// estimator may therefore be used, and in Stata that interface already
// exists: every estimation command posts e(b) and e(V).
// =====================================================================

struct phdid_data {
    real colvector tau
    real matrix Sigma
    real matrix Omega
    real colvector w
    real colvector g
    real colvector t
    string colvector lab
    real scalar K
    real scalar nunits
    real scalar hasgt
}

struct phdid_data scalar phdid_make(real colvector tau, real matrix Sigma,
                                    string colvector lab, real colvector w,
                                    real colvector g, real colvector t,
                                    real scalar nunits)
{
    struct phdid_data scalar d
    real scalar sw

    d.tau   = tau
    d.Sigma = (Sigma + Sigma') / 2
    d.Omega = phdid_spdinv(d.Sigma, "first-stage covariance")
    d.lab   = lab
    d.K     = rows(tau)
    d.nunits = nunits
    d.g = g
    d.t = t
    d.hasgt = (rows(g) == d.K & missing(g) == 0 & missing(t) == 0)

    if (rows(w) != d.K) {
        d.w = J(d.K, 1, 1 / d.K)
    }
    else {
        sw = colsum(w)
        if (sw <= 0) {
            d.w = J(d.K, 1, 1 / d.K)
        }
        else {
            d.w = w / sw
        }
    }
    return(d)
}

// Is the covariance diagonal?  Remark 1 of the paper: treating the
// cohort-time estimates as independent understates the posterior
// variance of linear aggregates and misstates co-clustering
// probabilities, so a diagonal Sigma is worth flagging.
real scalar phdid_isdiag(real matrix S)
{
    real matrix off
    real scalar sc

    off = S - diag(diagonal(S))
    sc  = max(abs(diagonal(S)))
    if (sc == 0) {
        return(1)
    }
    return((max(abs(off)) / sc) < 1e-10)
}

// ---------------------------------------------------------------------
// Persistent store.  R returns a ph_data object the user holds onto;
// Stata has no such thing, so the assembled data live in an external
// Mata global, in the same spirit as a posted estimation result.
// ---------------------------------------------------------------------

void phdid_setdata(struct phdid_data scalar d)
{
    external struct phdid_data scalar PHDID_D
    PHDID_D = d
}

struct phdid_data scalar phdid_getdata()
{
    external struct phdid_data scalar PHDID_D

    if (PHDID_D.K == . | PHDID_D.K == 0) {
        errprintf("no phdid data in memory\n")
        errprintf("{p 4 4 2}Run {bf:phdid setup} after a first-stage ")
        errprintf("estimation command, or supply b() and v() directly.{p_end}\n")
        exit(error(301))
    }
    return(PHDID_D)
}

// Assemble from Stata matrices and push into the store.  Returns the
// number of cells so the ado layer can report it.
real scalar phdid_setup(string scalar bname, string scalar vname,
                        string scalar labname, string scalar wname,
                        string scalar gname, string scalar tname,
                        real scalar nunits)
{
    struct phdid_data scalar d
    real colvector tau, w, g, t
    real matrix S
    string colvector lab
    real scalar K, i

    tau = st_matrix(bname)'
    S   = st_matrix(vname)
    K   = rows(tau)

    if (rows(S) != K | cols(S) != K) {
        errprintf("coefficient vector and covariance matrix do not conform\n")
        exit(error(503))
    }

    lab = J(K, 1, "")
    if (labname != "") {
        lab = tokens(st_global(labname))'
    }
    if (rows(lab) != K) {
        lab = J(K, 1, "")
        for (i = 1; i <= K; i++) {
            lab[i] = "cell" + strofreal(i)
        }
    }

    w = J(0, 1, .)
    if (wname != "") {
        w = st_matrix(wname)'
    }
    g = J(0, 1, .)
    t = J(0, 1, .)
    if (gname != "") {
        g = st_matrix(gname)'
        t = st_matrix(tname)'
    }

    d = phdid_make(tau, S, lab, w, g, t, nunits)
    phdid_setdata(d)
    return(K)
}

// Scalar accessors, so the ado layer never has to chain onto a
// struct-returning call (Mata rejects phdid_getdata().K).
real scalar phdid_K()
{
    struct phdid_data scalar d
    d = phdid_getdata()
    return(d.K)
}

real scalar phdid_hasgt()
{
    struct phdid_data scalar d
    d = phdid_getdata()
    return(d.hasgt)
}

// Report what is in the store.
void phdid_describe()
{
    struct phdid_data scalar d
    real matrix res
    real scalar k

    d = phdid_getdata()
    res = J(d.K, 3, .)
    for (k = 1; k <= d.K; k++) {
        res[k, 1] = d.tau[k]
        res[k, 2] = sqrt(d.Sigma[k, k])
        res[k, 3] = d.w[k]
    }
    st_matrix("r(cells)", res)
    st_matrixcolstripe("r(cells)", (J(3,1,""), ("estimate" \ "se" \ "weight")))
    st_matrixrowstripe("r(cells)", (J(d.K,1,""), d.lab))
    st_numscalar("r(K)", d.K)
    st_numscalar("r(isdiag)", phdid_isdiag(d.Sigma))
    st_numscalar("r(nunits)", d.nunits)
    st_numscalar("r(snr_desc)", sqrt(quadcross(d.tau :- mean(d.tau), d.tau :- mean(d.tau)) / (d.K - 1)) / mm_median_phdid(sqrt(diagonal(d.Sigma))))
    st_global("r(labels)", invtokens(d.lab'))
}

// =====================================================================
// L0 DRIVER
// =====================================================================

void phdid_l0_run(string scalar sel, real scalar lambda, real scalar mfix,
                  real scalar nbic)
{
    struct phdid_data scalar d
    struct phdid_path scalar P
    struct phdid_fit scalar f
    real colvector dev, bic, cgp, zsel
    real matrix res, pathm
    real scalar m, msel, bb, K, k, att, vat

    d = phdid_getdata()
    K = d.K
    P = phdid_agglom(d.tau, d.Omega, .)

    dev = J(K, 1, .)
    bic = J(K, 1, .)
    cgp = J(K, 1, .)
    for (m = 1; m <= K; m++) {
        f = phdid_gls(d.tau, d.Omega, P.parts[., m])
        dev[m] = f.dev
        bic[m] = phdid_bic(f.dev, m, nbic)
        cgp[m] = phdid_crosspairs(P.parts[., m])
    }

    if (sel == "bic") {
        msel = .
        bb = .
        for (m = 1; m <= K; m++) {
            if (bic[m] < bb) {
                bb = bic[m]
                msel = m
            }
        }
    }
    else if (sel == "lambda") {
        msel = phdid_select_lambda(P, lambda, K)
    }
    else {
        msel = mfix
    }

    zsel = P.parts[., msel]
    f = phdid_gls(d.tau, d.Omega, zsel)

    // per-cell results: flexible, se, group, grouped, se, variance ratio
    res = J(K, 6, .)
    for (k = 1; k <= K; k++) {
        res[k, 1] = d.tau[k]
        res[k, 2] = sqrt(d.Sigma[k, k])
        res[k, 3] = f.z[k]
        res[k, 4] = f.tau[k]
        res[k, 5] = sqrt(f.V[k, k])
        res[k, 6] = f.V[k, k] / d.Sigma[k, k]
    }
    pathm = (1::K), dev, bic, cgp, P.mergecost, P.pairprod

    att = quadcross(d.w, f.tau)
    vat = quadcross(d.w, f.V * d.w)

    st_matrix("r(cells)", res)
    st_matrixcolstripe("r(cells)", (J(6,1,""),
        ("flexible" \ "se_flex" \ "group" \ "grouped" \ "se" \ "var_ratio")))
    st_matrixrowstripe("r(cells)", (J(K,1,""), d.lab))
    st_matrix("r(path)", pathm)
    st_matrixcolstripe("r(path)", (J(6,1,""),
        ("m" \ "deviance" \ "bic" \ "crosspairs" \ "mergecost" \ "pairprod")))
    st_matrix("r(b)", f.tau')
    st_matrixcolstripe("r(b)", (J(K,1,""), d.lab))
    st_matrix("r(V)", f.V)
    st_matrixcolstripe("r(V)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(V)", (J(K,1,""), d.lab))
    st_matrix("r(phi)", f.phi')
    st_matrix("r(partition)", f.z')
    st_numscalar("r(m)", f.m)
    st_numscalar("r(deviance)", f.dev)
    st_numscalar("r(bic)", bic[msel])
    st_numscalar("r(att)", att)
    st_numscalar("r(att_se)", sqrt(vat))
    st_numscalar("r(K)", K)
    st_numscalar("r(nbic)", nbic)
    st_global("r(partition_str)", phdid_fmtpart(f.z, d.lab))
}

// =====================================================================
// SPECIFICATION TEST  (paper Sections 5.1 and 5.2.3)
//
// Provenance matters here, because this package ships alongside the
// paper.  Only the placebo gauge and the within-cell randomization test
// are the paper's procedure.  The common-effect chi-squared test and
// the covariance-based dispersion decomposition were added for the
// software, because the placebo route needs pre-treatment estimates
// that many designs -- including the paper's own first application --
// cannot supply.
// =====================================================================

void phdid_test_run(real colvector pre, real scalar nsim, real scalar center)
{
    struct phdid_data scalar d
    struct phdid_fit scalar pl
    real matrix M
    real colvector postc, prec, drawn, ref
    real scalar K, obsdisp, expdisp, hetvar, share, snr, snrdesc
    real scalar spost, spre, hv2, sh2, obs, b, i, pval, nge

    d = phdid_getdata()
    K = d.K

    // 1. common-effect test on the pooled GLS deviance ~ chi2(K-1)
    pl = phdid_gls(d.tau, d.Omega, J(K, 1, 1))
    st_numscalar("r(chi2)", pl.dev)
    st_numscalar("r(df)", K - 1)
    st_numscalar("r(p)", chi2tail(K - 1, pl.dev))

    // 2. excess dispersion against what the covariance predicts
    M = I(K) :- J(K, K, 1 / K)
    obsdisp = quadcross(d.tau :- mean(d.tau), d.tau :- mean(d.tau)) / (K - 1)
    expdisp = trace(M * d.Sigma * M) / (K - 1)
    hetvar  = max((0, obsdisp - expdisp))
    if (obsdisp > 0) {
        share = hetvar / obsdisp
    }
    else {
        share = 0
    }
    snr = sqrt(obsdisp) / sqrt(expdisp)
    snrdesc = sqrt(obsdisp) / mm_median_phdid(sqrt(diagonal(d.Sigma)))

    st_numscalar("r(disp_obs)", obsdisp)
    st_numscalar("r(disp_exp)", expdisp)
    st_numscalar("r(het_var)", hetvar)
    st_numscalar("r(het_share)", share)
    st_numscalar("r(snr)", snr)
    st_numscalar("r(snr_desc)", snrdesc)
    st_numscalar("r(has_placebo)", 0)

    // 3. placebo gauge and within-cell randomization test
    if (rows(pre) == K) {
        spost = sqrt(quadcross(d.tau :- mean(d.tau), d.tau :- mean(d.tau)) / (K - 1))
        spre  = sqrt(quadcross(pre :- mean(pre), pre :- mean(pre)) / (K - 1))
        hv2 = max((0, spost^2 - spre^2))
        if (spost > 0) {
            sh2 = hv2 / spost^2
        }
        else {
            sh2 = 0
        }
        if (center == 1) {
            postc = d.tau :- mean(d.tau)
            prec  = pre :- mean(pre)
        }
        else {
            postc = d.tau
            prec  = pre
        }
        obs = sqrt(quadcross(postc :- mean(postc), postc :- mean(postc)) / (K - 1))
        ref = J(nsim, 1, .)
        for (b = 1; b <= nsim; b++) {
            drawn = J(K, 1, .)
            for (i = 1; i <= K; i++) {
                if (runiform(1, 1) < 0.5) {
                    drawn[i] = prec[i]
                }
                else {
                    drawn[i] = postc[i]
                }
            }
            ref[b] = sqrt(quadcross(drawn :- mean(drawn), drawn :- mean(drawn)) / (K - 1))
        }
        nge = colsum(ref :>= obs)
        pval = (1 + nge) / (nsim + 1)

        st_numscalar("r(has_placebo)", 1)
        st_numscalar("r(s_pre)", spre)
        st_numscalar("r(s_post)", spost)
        st_numscalar("r(pl_snr)", spost / spre)
        st_numscalar("r(pl_het_share)", sh2)
        st_numscalar("r(pl_p)", pval)
        st_numscalar("r(pl_obs)", obs)
        st_numscalar("r(pl_refmed)", mm_median_phdid(ref))
        st_numscalar("r(nsim)", nsim)
        st_numscalar("r(center)", center)
    }
}

// Median without depending on moremata, which is not a Stata Journal
// dependency we want to take on.
real scalar mm_median_phdid(real colvector x)
{
    real colvector s
    real scalar n

    s = sort(x, 1)
    n = rows(s)
    if (mod(n, 2) == 1) {
        return(s[(n + 1) / 2])
    }
    return((s[n / 2] + s[n / 2 + 1]) / 2)
}

// =====================================================================
// FIT UNDER A PARTITION THE USER SUPPLIES
//
// The estimator's two halves are separable: Proposition 1 gives the
// restricted GLS fit for any partition, and the search in Appendix B is
// only one way to choose one.  A partition that comes from theory, from
// an earlier sample, or from the posterior can be fitted directly here,
// and its deviance compared with the selected one.
// =====================================================================

void phdid_fit_run(real colvector zin)
{
    struct phdid_data scalar d
    struct phdid_fit scalar f
    real matrix res
    real scalar K, k, att, vat

    d = phdid_getdata()
    K = d.K
    if (rows(zin) != K) {
        errprintf("partition() has %f entries but there are %f cells\n",
                  rows(zin), K)
        exit(error(503))
    }
    if (min(zin) < 1 | missing(zin) > 0) {
        errprintf("partition() must give every cell a positive group number\n")
        exit(error(198))
    }

    f = phdid_gls(d.tau, d.Omega, zin)

    res = J(K, 6, .)
    for (k = 1; k <= K; k++) {
        res[k, 1] = d.tau[k]
        res[k, 2] = sqrt(d.Sigma[k, k])
        res[k, 3] = f.z[k]
        res[k, 4] = f.tau[k]
        res[k, 5] = sqrt(f.V[k, k])
        res[k, 6] = f.V[k, k] / d.Sigma[k, k]
    }
    att = quadcross(d.w, f.tau)
    vat = quadcross(d.w, f.V * d.w)

    st_matrix("r(cells)", res)
    st_matrixcolstripe("r(cells)", (J(6,1,""),
        ("flexible" \ "se_flex" \ "group" \ "grouped" \ "se" \ "var_ratio")))
    st_matrixrowstripe("r(cells)", (J(K,1,""), d.lab))
    st_matrix("r(b)", f.tau')
    st_matrixcolstripe("r(b)", (J(K,1,""), d.lab))
    st_matrix("r(V)", f.V)
    st_matrixcolstripe("r(V)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(V)", (J(K,1,""), d.lab))
    st_matrix("r(phi)", f.phi')
    st_matrix("r(partition)", f.z')
    st_matrixcolstripe("r(partition)", (J(K,1,""), d.lab))
    st_numscalar("r(m)", f.m)
    st_numscalar("r(deviance)", f.dev)
    st_numscalar("r(df)", K - f.m)
    st_numscalar("r(p)", chi2tail(K - f.m, f.dev))
    st_numscalar("r(att)", att)
    st_numscalar("r(att_se)", sqrt(vat))
    st_numscalar("r(K)", K)
    st_global("r(partition_str)", phdid_fmtpart(f.z, d.lab))
}

end
