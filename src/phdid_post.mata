*! version 0.1.0  phdid post-estimation core
*! Rohan Wagle and Parush Arora
// =====================================================================
// Post-estimation: aggregation, sampler diagnostics, point partitions,
// the exact-vs-diagonal covariance check and the two regularisation
// paths.  Everything here reads the draws left in memory by
// phdid_gibbs(), or the posted e(b) and e(V) of a conditional fit.
// =====================================================================

version 17.0
mata: mata set matastrict on
mata:

// Percentile interval of a draw vector at a confidence level in percent.
real rowvector phdid_pctl(real colvector x, real scalar level)
{
    real colvector q
    real scalar n, a

    n = rows(x)
    a = (1 - level / 100) / 2
    q = sort(x, 1)
    return((q[max((1, ceil(a * n)))], q[min((n, ceil((1 - a) * n)))]))
}

real scalar phdid_sd(real colvector x)
{
    real scalar n

    n = rows(x)
    if (n < 2) {
        return(.)
    }
    return(sqrt(quadcross(x :- mean(x), x :- mean(x)) / (n - 1)))
}

// =====================================================================
// AGGREGATION
//
// Every estimand the package reports is a fixed linear map of the cell
// effects, so an aggregation is a matrix A with one row per term.  For
// the Bayes fit the map is applied to each posterior draw, which is what
// makes the reported interval include uncertainty about the partition
// itself; for the l0 fit it is applied to the posted e(b) and e(V),
// which conditions on the selected partition being right.
// =====================================================================

struct phdid_agg {
    real matrix A
    string colvector lab
}

struct phdid_agg scalar phdid_aggmat(struct phdid_data scalar d,
                                     string scalar type)
{
    struct phdid_agg scalar G
    real colvector key, u, wk
    real scalar n, i, sw
    string scalar pre

    if (type == "overall") {
        G.A   = d.w'
        G.lab = "overall"
        return(G)
    }

    if (d.hasgt == 0) {
        errprintf("cohort and period identifiers are not available\n")
        errprintf("{p 4 4 2}The %s aggregation needs to know which cohort ", type)
        errprintf("and period each cell belongs to.  {bf:phdid setup} ")
        errprintf("recovers these from the first-stage coefficient names; ")
        errprintf("supply {bf:cells()} when they cannot be parsed.{p_end}\n")
        exit(error(198))
    }

    pre = ""
    if (type == "dynamic") {
        key = d.t - d.g
        pre = "e="
    }
    else if (type == "group") {
        key = d.g
        pre = "g="
    }
    else if (type == "calendar") {
        key = d.t
        pre = "t="
    }
    else {
        errprintf("unknown aggregation %s\n", type)
        exit(error(198))
    }

    u = uniqrows(key)
    n = rows(u)
    G.A   = J(n, d.K, 0)
    G.lab = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        wk = d.w :* (key :== u[i])
        sw = colsum(wk)
        if (sw > 0) {
            G.A[i, .] = (wk / sw)'
        }
        G.lab[i] = pre + strofreal(u[i])
    }
    return(G)
}

void phdid_agg_run(string scalar type, real scalar level, string scalar src)
{
    external real matrix PHDID_TAUDRAWS

    struct phdid_data scalar d
    struct phdid_agg scalar G
    real matrix AD, res, Vout, V
    real colvector b, est, se
    real rowvector ci
    real scalar n, i, nd, zq

    d = phdid_getdata()
    G = phdid_aggmat(d, type)
    n = rows(G.A)
    res = J(n, 4, .)

    if (src == "draws") {
        if (rows(PHDID_TAUDRAWS) == 0) {
            errprintf("no posterior draws in memory\n")
            errprintf("{p 4 4 2}Run {bf:phdid bayes} first.{p_end}\n")
            exit(error(301))
        }
        AD = PHDID_TAUDRAWS * G.A'
        nd = rows(AD)
        for (i = 1; i <= n; i++) {
            res[i, 1] = mean(AD[., i])
            res[i, 2] = phdid_sd(AD[., i])
            ci = phdid_pctl(AD[., i], level)
            res[i, 3] = ci[1]
            res[i, 4] = ci[2]
        }
        Vout = variance(AD)
        st_numscalar("r(ndraws)", nd)
    }
    else {
        b = st_matrix("e(b)")'
        V = st_matrix("e(V)")
        if (rows(b) != d.K) {
            errprintf("the stored fit has %f coefficients but the data in ", rows(b))
            errprintf("memory have %f cells\n", d.K)
            exit(error(503))
        }
        est  = G.A * b
        Vout = G.A * V * G.A'
        se   = sqrt(diagonal(Vout))
        zq   = invnormal(1 - (1 - level / 100) / 2)
        for (i = 1; i <= n; i++) {
            res[i, 1] = est[i]
            res[i, 2] = se[i]
            res[i, 3] = est[i] - zq * se[i]
            res[i, 4] = est[i] + zq * se[i]
        }
    }

    st_matrix("r(table)", res)
    st_matrixcolstripe("r(table)", (J(4,1,""),
        ("estimate" \ "se" \ "lower" \ "upper")))
    st_matrixrowstripe("r(table)", (J(n,1,""), G.lab))
    st_matrix("r(b)", res[., 1]')
    st_matrixcolstripe("r(b)", (J(n,1,""), G.lab))
    st_matrix("r(V)", Vout)
    st_matrixcolstripe("r(V)", (J(n,1,""), G.lab))
    st_matrixrowstripe("r(V)", (J(n,1,""), G.lab))
    st_matrix("r(A)", G.A)
    st_matrixrowstripe("r(A)", (J(n,1,""), G.lab))
    st_matrixcolstripe("r(A)", (J(d.K,1,""), d.lab))
    st_numscalar("r(nterms)", n)
    st_global("r(aggregation)", type)
}

// =====================================================================
// SAMPLER DIAGNOSTICS
//
// Split-chain Gelman-Rubin on each cell effect and on the number of
// groups, with effective sample sizes from Geyer's initial positive
// sequence.  The partition label itself is not a continuous parameter,
// so the number of groups is the summary that carries the mixing: a
// chain stuck in one corner of partition space shows up there.
// =====================================================================

// Autocovariance of x at lag t, divided by n.
real scalar phdid_acov(real colvector x, real scalar t, real scalar xbar)
{
    real scalar n

    n = rows(x)
    if (t >= n) {
        return(0)
    }
    return(quadcross(x[|1 \ n - t|] :- xbar, x[|1 + t \ n|] :- xbar) / n)
}

// Split-Rhat and effective sample size for one column of draws.
real rowvector phdid_rhat1(real colvector x, real colvector chain,
                           real scalar maxlag)
{
    real matrix seq
    real colvector ids, mn, vr, xm
    real scalar nch, nseq, half, i, s, W, B, vhat, rhat, t, rho, tsum, ess
    real scalar acv, prev, pair, lmax

    ids  = uniqrows(chain)
    nch  = rows(ids)
    half = floor(rows(x) / nch / 2)
    if (half < 4) {
        return((., ., .))
    }

    nseq = 2 * nch
    seq  = J(half, nseq, .)
    s    = 0
    for (i = 1; i <= nch; i++) {
        xm = select(x, chain :== ids[i])
        s  = s + 1
        seq[., s] = xm[|1 \ half|]
        s  = s + 1
        seq[., s] = xm[|rows(xm) - half + 1 \ rows(xm)|]
    }

    mn = J(nseq, 1, .)
    vr = J(nseq, 1, .)
    for (s = 1; s <= nseq; s++) {
        mn[s] = mean(seq[., s])
        vr[s] = phdid_sd(seq[., s])^2
    }
    W = mean(vr)
    B = phdid_sd(mn)^2
    if (W <= 0 | W >= .) {
        return((., ., half * nseq))
    }
    vhat = ((half - 1) / half) * W + B
    rhat = sqrt(vhat / W)

    // Geyer's initial positive sequence, truncated at the first negative
    // sum of an adjacent pair of autocorrelations.
    lmax = min((maxlag, half - 2))
    tsum = 0
    prev = .
    for (t = 1; t <= lmax; t++) {
        acv = 0
        for (s = 1; s <= nseq; s++) {
            acv = acv + phdid_acov(seq[., s], t, mn[s])
        }
        acv = acv / nseq
        rho = 1 - (W - acv) / vhat
        if (mod(t, 2) == 1) {
            prev = rho
        }
        else {
            pair = prev + rho
            if (pair < 0) {
                break
            }
            tsum = tsum + pair
        }
    }
    ess = (half * nseq) / (1 + 2 * tsum)
    if (ess > half * nseq) {
        ess = half * nseq
    }
    return((rhat, ess, half * nseq))
}

void phdid_rhat_run(real scalar maxlag)
{
    external real matrix PHDID_TAUDRAWS
    external real colvector PHDID_MDRAWS
    external real colvector PHDID_CHAIN

    struct phdid_data scalar d
    real matrix res
    real rowvector r1
    real scalar K, k, nch

    d = phdid_getdata()
    K = d.K
    if (rows(PHDID_TAUDRAWS) == 0) {
        errprintf("no posterior draws in memory\n")
        errprintf("{p 4 4 2}Run {bf:phdid bayes} first.{p_end}\n")
        exit(error(301))
    }
    nch = rows(uniqrows(PHDID_CHAIN))

    res = J(K + 1, 3, .)
    for (k = 1; k <= K; k++) {
        r1 = phdid_rhat1(PHDID_TAUDRAWS[., k], PHDID_CHAIN, maxlag)
        res[k, .] = r1
    }
    res[K + 1, .] = phdid_rhat1(PHDID_MDRAWS, PHDID_CHAIN, maxlag)

    st_matrix("r(table)", res)
    st_matrixcolstripe("r(table)", (J(3,1,""), ("rhat" \ "ess" \ "n_split")))
    st_matrixrowstripe("r(table)", (J(K + 1,1,""), (d.lab \ "groups (m)")))
    st_numscalar("r(nchains)", nch)
    st_numscalar("r(rhat_max)", max(res[., 1]))
    st_numscalar("r(ess_min)", min(res[., 2]))
}

// =====================================================================
// POINT PARTITION FROM THE POSTERIOR
//
// The posterior mode over partitions is a poor summary -- with K cells
// no single partition holds much mass -- so the package minimises an
// expected loss over the sampled partitions instead: Binder's loss,
// which is the L1 distance to the co-clustering matrix, or the
// variation of information of Wade and Ghahramani (2018), which
// respects the lattice structure of partition space.
// =====================================================================

real scalar phdid_vi(real colvector z1, real colvector z2)
{
    real matrix nab
    real colvector na, nb
    real scalar K, m1, m2, a, b, vi, pa, pb, pab

    K  = rows(z1)
    m1 = max(z1)
    m2 = max(z2)
    nab = J(m1, m2, 0)
    for (a = 1; a <= K; a++) {
        nab[z1[a], z2[a]] = nab[z1[a], z2[a]] + 1
    }
    na = rowsum(nab)
    nb = colsum(nab)'

    vi = 0
    for (a = 1; a <= m1; a++) {
        pa = na[a] / K
        if (pa > 0) {
            vi = vi - pa * ln(pa)
        }
    }
    for (b = 1; b <= m2; b++) {
        pb = nb[b] / K
        if (pb > 0) {
            vi = vi - pb * ln(pb)
        }
    }
    for (a = 1; a <= m1; a++) {
        for (b = 1; b <= m2; b++) {
            if (nab[a, b] > 0) {
                pab = nab[a, b] / K
                pa  = na[a] / K
                pb  = nb[b] / K
                vi  = vi - 2 * pab * ln(pab / (pa * pb))
            }
        }
    }
    return(vi)
}

void phdid_point_run(string scalar loss)
{
    external real matrix PHDID_ZDRAWS

    struct phdid_data scalar d
    struct phdid_fit scalar f
    real matrix Zu, cc, Zr
    real colvector cnt, obj, zc, zs, att
    real scalar K, nd, nu, i, j, k, L, tot, best, bi

    d = phdid_getdata()
    K = d.K
    if (rows(PHDID_ZDRAWS) == 0) {
        errprintf("no posterior draws in memory\n")
        errprintf("{p 4 4 2}Run {bf:phdid bayes} first.{p_end}\n")
        exit(error(301))
    }

    nd = rows(PHDID_ZDRAWS)
    Zr = J(nd, K, .)
    for (i = 1; i <= nd; i++) {
        Zr[i, .] = phdid_relabel(PHDID_ZDRAWS[i, .]')'
    }
    Zu  = uniqrows(Zr)
    nu  = rows(Zu)
    cnt = J(nu, 1, 0)
    for (i = 1; i <= nd; i++) {
        for (j = 1; j <= nu; j++) {
            if (Zr[i, .] == Zu[j, .]) {
                cnt[j] = cnt[j] + 1
                break
            }
        }
    }

    cc = J(K, K, 0)
    for (i = 1; i <= nd; i++) {
        for (j = 1; j <= K; j++) {
            for (k = 1; k <= K; k++) {
                if (Zr[i, j] == Zr[i, k]) {
                    cc[j, k] = cc[j, k] + 1
                }
            }
        }
    }
    cc = cc / nd

    obj = J(nu, 1, .)
    for (i = 1; i <= nu; i++) {
        zc = Zu[i, .]'
        if (loss == "binder") {
            tot = 0
            for (j = 1; j <= K - 1; j++) {
                for (k = j + 1; k <= K; k++) {
                    tot = tot + abs((zc[j] == zc[k]) - cc[j, k])
                }
            }
            obj[i] = tot
        }
        else {
            tot = 0
            for (j = 1; j <= nu; j++) {
                zs  = Zu[j, .]'
                tot = tot + cnt[j] * phdid_vi(zc, zs)
            }
            obj[i] = tot / nd
        }
    }

    best = .
    bi   = .
    for (i = 1; i <= nu; i++) {
        if (obj[i] < best) {
            best = obj[i]
            bi   = i
        }
    }

    zc = Zu[bi, .]'
    f  = phdid_gls(d.tau, d.Omega, zc)
    att = J(1, 1, quadcross(d.w, f.tau))

    st_matrix("r(partition)", zc')
    st_matrixcolstripe("r(partition)", (J(K,1,""), d.lab))
    st_matrix("r(b)", f.tau')
    st_matrixcolstripe("r(b)", (J(K,1,""), d.lab))
    st_matrix("r(V)", f.V)
    st_matrixcolstripe("r(V)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(V)", (J(K,1,""), d.lab))
    st_matrix("r(coclust)", cc)
    st_matrixcolstripe("r(coclust)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(coclust)", (J(K,1,""), d.lab))
    st_numscalar("r(m)", f.m)
    st_numscalar("r(loss_value)", best)
    st_numscalar("r(n_unique)", nu)
    st_numscalar("r(post_prob)", cnt[bi] / nd)
    st_numscalar("r(att)", att[1, 1])
    st_global("r(loss)", loss)
    st_global("r(partition_str)", phdid_fmtpart(zc, d.lab))
}

// =====================================================================
// EXACT VERSUS DIAGONAL COVARIANCE  (paper Remark 1)
//
// Runs the same sampler twice, once with the exact collapsed marginal
// likelihood and once with the diagonal shortcut, and reports how much
// the shortcut moves.  Cells share units and share the two-way fixed
// effects, so the shortcut is not a harmless speed-up: it understates
// the posterior variance of aggregates and misstates individual
// co-clustering probabilities.
// =====================================================================

void phdid_covcheck_run(real scalar alpha, real scalar mu0in,
                        real scalar s0sqin, real scalar iters,
                        real scalar burn, real scalar thin,
                        real scalar chains, string scalar init)
{
    struct phdid_data scalar d
    struct phdid_draws scalar De
    struct phdid_draws scalar Dd
    struct phdid_agg scalar G
    real rowvector base
    real matrix cce, ccd, res
    real colvector ae, ad
    real scalar K, mu0, s0sq, i, j, mx, sde, sdd

    d = phdid_getdata()
    K = d.K
    base = phdid_basemeasure(d, mu0in, s0sqin)
    mu0  = base[1]
    s0sq = base[2]

    De = phdid_gibbs(d, alpha, mu0, s0sq, iters, burn, thin, chains, 1, init)
    Dd = phdid_gibbs(d, alpha, mu0, s0sq, iters, burn, thin, chains, 0, init)

    cce = J(K, K, 0)
    ccd = J(K, K, 0)
    for (i = 1; i <= De.n; i++) {
        for (j = 1; j <= K; j++) {
            cce[j, .] = cce[j, .] + (De.z[i, j] :== De.z[i, .])
        }
    }
    for (i = 1; i <= Dd.n; i++) {
        for (j = 1; j <= K; j++) {
            ccd[j, .] = ccd[j, .] + (Dd.z[i, j] :== Dd.z[i, .])
        }
    }
    cce = cce / De.n
    ccd = ccd / Dd.n

    G  = phdid_aggmat(d, "overall")
    ae = De.tau * G.A'
    ad = Dd.tau * G.A'
    sde = phdid_sd(ae)
    sdd = phdid_sd(ad)

    mx = 0
    for (i = 1; i <= K; i++) {
        for (j = 1; j <= K; j++) {
            if (abs(cce[i, j] - ccd[i, j]) > mx) {
                mx = abs(cce[i, j] - ccd[i, j])
            }
        }
    }

    res = J(K, 2, .)
    for (i = 1; i <= K; i++) {
        res[i, 1] = mean(De.tau[., i])
        res[i, 2] = mean(Dd.tau[., i])
    }

    st_matrix("r(cells)", res)
    st_matrixcolstripe("r(cells)", (J(2,1,""), ("exact" \ "diagonal")))
    st_matrixrowstripe("r(cells)", (J(K,1,""), d.lab))
    st_matrix("r(cc_exact)", cce)
    st_matrix("r(cc_diag)", ccd)
    st_numscalar("r(cc_maxdiff)", mx)
    st_numscalar("r(m_exact)", mean(De.m))
    st_numscalar("r(m_diag)", mean(Dd.m))
    st_numscalar("r(att_exact)", mean(ae))
    st_numscalar("r(att_diag)", mean(ad))
    st_numscalar("r(sd_exact)", sde)
    st_numscalar("r(sd_diag)", sdd)
    st_numscalar("r(var_ratio)", (sde / sdd)^2)
    st_numscalar("r(isdiag)", phdid_isdiag(d.Sigma))
}

// =====================================================================
// REGULARISATION PATHS
//
// How much of what is reported is the prior talking?  Both paths answer
// that question by re-estimating across the tuning parameter.  The
// per-cell view is usually the informative one: the overall ATT is
// robust to over-pooling, so a flat aggregate path can hide cells that
// move a great deal.
// =====================================================================

void phdid_alphapath_run(real rowvector alphas, real scalar mu0in,
                         real scalar s0sqin, real scalar iters,
                         real scalar burn, real scalar thin,
                         real scalar chains, real scalar exact,
                         string scalar init, real scalar level,
                         string scalar type)
{
    struct phdid_data scalar d
    struct phdid_draws scalar D
    struct phdid_agg scalar G
    real rowvector base, ci
    real matrix res
    real colvector agg
    real scalar K, mu0, s0sq, na, i, k

    d = phdid_getdata()
    K = d.K
    base = phdid_basemeasure(d, mu0in, s0sqin)
    mu0  = base[1]
    s0sq = base[2]
    G    = phdid_aggmat(d, "overall")
    na   = cols(alphas)

    if (type == "cells") {
        res = J(na, 2 + K, .)
    }
    else {
        res = J(na, 6, .)
    }

    for (i = 1; i <= na; i++) {
        D = phdid_gibbs(d, alphas[i], mu0, s0sq, iters, burn, thin, chains,
                        exact, init)
        agg = D.tau * G.A'
        res[i, 1] = alphas[i]
        res[i, 2] = mean(D.m)
        if (type == "cells") {
            for (k = 1; k <= K; k++) {
                res[i, 2 + k] = mean(D.tau[., k])
            }
        }
        else {
            ci = phdid_pctl(agg, level)
            res[i, 3] = mean(agg)
            res[i, 4] = phdid_sd(agg)
            res[i, 5] = ci[1]
            res[i, 6] = ci[2]
        }
    }

    st_matrix("r(path)", res)
    if (type == "cells") {
        st_matrixcolstripe("r(path)", (J(2 + K,1,""),
            ("alpha" \ "m_mean" \ d.lab)))
    }
    else {
        st_matrixcolstripe("r(path)", (J(6,1,""),
            ("alpha" \ "m_mean" \ "att" \ "sd" \ "lower" \ "upper")))
    }
    st_numscalar("r(npoints)", na)
    st_global("r(parameter)", "alpha")
    st_global("r(path_type)", type)
}

void phdid_lambdapath_run(real rowvector lams, real scalar level,
                          string scalar type, string scalar by)
{
    struct phdid_data scalar d
    struct phdid_path scalar P
    struct phdid_fit scalar f
    real matrix res
    real colvector zsel
    real scalar K, nl, i, k, m, att, vat, zq, lam

    d = phdid_getdata()
    K = d.K
    P = phdid_agglom(d.tau, d.Omega, .)
    zq = invnormal(1 - (1 - level / 100) / 2)

    if (by == "m") {
        nl = K
    }
    else {
        nl = cols(lams)
    }

    if (type == "cells") {
        res = J(nl, 2 + K, .)
    }
    else {
        res = J(nl, 6, .)
    }

    for (i = 1; i <= nl; i++) {
        if (by == "m") {
            m   = K - i + 1
            lam = P.mergecost[m] / P.pairprod[m]
        }
        else {
            lam = lams[i]
            m   = phdid_select_lambda(P, lam, K)
        }
        zsel = P.parts[., m]
        f    = phdid_gls(d.tau, d.Omega, zsel)
        att  = quadcross(d.w, f.tau)
        vat  = quadcross(d.w, f.V * d.w)

        res[i, 1] = lam
        res[i, 2] = f.m
        if (type == "cells") {
            for (k = 1; k <= K; k++) {
                res[i, 2 + k] = f.tau[k]
            }
        }
        else {
            res[i, 3] = att
            res[i, 4] = sqrt(vat)
            res[i, 5] = att - zq * sqrt(vat)
            res[i, 6] = att + zq * sqrt(vat)
        }
    }

    st_matrix("r(path)", res)
    if (type == "cells") {
        st_matrixcolstripe("r(path)", (J(2 + K,1,""),
            ("lambda" \ "m" \ d.lab)))
    }
    else {
        st_matrixcolstripe("r(path)", (J(6,1,""),
            ("lambda" \ "m" \ "att" \ "se" \ "lower" \ "upper")))
    }
    st_matrix("r(thresholds)", (P.mergecost :/ P.pairprod)')
    st_numscalar("r(npoints)", nl)
    st_global("r(parameter)", "lambda")
    st_global("r(path_type)", type)
    st_global("r(by)", by)
}

end
