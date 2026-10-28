*! phdid_bayes.mata -- Dirichlet process estimator for the phdid package
*! Rohan Wagle and Parush Arora
version 17.0

mata: mata set matastrict on

mata:

// =====================================================================
// COLLAPSED GIBBS SAMPLER  (paper Section 3; Neal 2000, Algorithm 3)
//
// Everything works in the two-stage coordinates of Lemma 1: the data
// enter only through the group-summed precision S = R'Omega R and the
// group-summed score u = R'Omega tau.
//
// The collapsed marginal likelihood of a partition is eq. (30), in the
// general-covariance form of eq. (54):
//
//   log p(y|P) = -1/2 log|I_m + s0sq*S| + s0sq/2 * u0'(I_m+s0sq*S)^-1 u0
//
// with u0 = u - mu0*S*1_m.  Terms constant across partitions drop,
// which is legitimate because R*1_m = 1_K for every partition, so the
// centred outcome does not depend on P.
// =====================================================================

real scalar phdid_logml(real colvector z, real matrix Omega,
                        real colvector Otau, real scalar mu0,
                        real scalar s0sq)
{
    real matrix R, S, A, L
    real colvector u, u0, sol
    real scalar m, ldet, qform

    R = phdid_indicator(z)
    S = quadcross(R, Omega * R)
    u = quadcross(R, Otau)
    m = cols(R)

    u0 = u - mu0 * (S * J(m, 1, 1))
    A  = I(m) + s0sq * S
    A  = (A + A') / 2

    L = cholesky(A)
    if (missing(L) > 0) {
        return(-9e300)
    }
    ldet = 2 * colsum(ln(diagonal(L)))
    sol  = cholsolve(A, u0)
    qform = quadcross(u0, sol)

    return(-0.5 * ldet + 0.5 * s0sq * qform)
}

// Group-effect posterior, eq. (31)/(55):
//   Sigma* = (S + I_m/s0sq)^-1 ,  mu* = Sigma*(u + mu0*1_m/s0sq)
// This is the ridge-regularised GLS estimator of the group effects.
real colvector phdid_draw_phi(real colvector z, real matrix Omega,
                              real colvector Otau, real scalar mu0,
                              real scalar s0sq)
{
    real matrix R, S, P, L
    real colvector u, mu, eps

    R = phdid_indicator(z)
    S = quadcross(R, Omega * R)
    u = quadcross(R, Otau)
    P = S + I(cols(R)) / s0sq
    P = (P + P') / 2

    mu = cholsolve(P, u :+ mu0 / s0sq)
    L  = cholesky(P)
    eps = rnormal(cols(R), 1, 0, 1)

    // x = mu + L'^-1 eps has covariance P^-1 when P = L L'
    return(mu + lusolve(L', eps))
}

// One sweep of cluster assignments.  For each cell: remove it, score
// every surviving cluster and a fresh one, draw.
real colvector phdid_sweep(real colvector zin, real matrix Omega,
                           real colvector Otau, real colvector tau,
                           real colvector prec, real scalar alpha,
                           real scalar mu0, real scalar s0sq,
                           real scalar exact)
{
    real colvector z, zminus, labs, logp, pr, trial, members
    real scalar K, k, j, ncl, nmem, pick, mx, tot, cum, uu, Pc, mc, i

    z = zin
    K = rows(z)
    for (k = 1; k <= K; k++) {
        zminus = z
        zminus[k] = 0
        labs = uniqrows(select(zminus, zminus :> 0))
        ncl  = rows(labs)
        logp = J(ncl + 1, 1, .)

        for (j = 1; j <= ncl; j++) {
            members = select((1::K), zminus :== labs[j])
            nmem = rows(members)
            if (exact == 1) {
                trial = zminus
                trial[k] = labs[j]
                logp[j] = ln(nmem) + phdid_logml(phdid_relabel(trial), Omega, Otau, mu0, s0sq)
            }
            else {
                Pc = 1 / s0sq
                mc = mu0 / s0sq
                for (i = 1; i <= nmem; i++) {
                    Pc = Pc + prec[members[i]]
                    mc = mc + prec[members[i]] * tau[members[i]]
                }
                mc = mc / Pc
                logp[j] = ln(nmem) + lnnormalden(tau[k], mc, sqrt(1 / prec[k] + 1 / Pc))
            }
        }

        if (exact == 1) {
            trial = zminus
            trial[k] = max(labs \ 0) + 1
            logp[ncl + 1] = ln(alpha) + phdid_logml(phdid_relabel(trial), Omega, Otau, mu0, s0sq)
        }
        else {
            logp[ncl + 1] = ln(alpha) + lnnormalden(tau[k], mu0, sqrt(1 / prec[k] + s0sq))
        }

        mx = max(logp)
        pr = exp(logp :- mx)
        tot = colsum(pr)
        if (tot <= 0) {
            continue
        }
        pr = pr / tot

        uu = runiform(1, 1)
        cum = 0
        pick = ncl + 1
        for (j = 1; j <= ncl + 1; j++) {
            cum = cum + pr[j]
            if (uu <= cum) {
                pick = j
                break
            }
        }

        z = zminus
        if (pick <= ncl) {
            z[k] = labs[pick]
        }
        else {
            z[k] = max(labs \ 0) + 1
        }
        z = phdid_relabel(z)
    }
    return(z)
}

// Dispersed starting partitions, as in the paper's Appendix E: all
// singletons, one pooled group, and two random partitions.
real colvector phdid_init(string scalar how, real scalar K, real scalar chain)
{
    real scalar wch

    if (how == "singletons") {
        return((1::K))
    }
    if (how == "pooled") {
        return(J(K, 1, 1))
    }
    if (how == "random") {
        return(phdid_relabel(ceil(runiform(K, 1) * max((2, floor(K / 2))))))
    }
    wch = mod(chain - 1, 4) + 1
    if (wch == 1) {
        return((1::K))
    }
    if (wch == 2) {
        return(J(K, 1, 1))
    }
    if (wch == 3) {
        return(phdid_relabel(ceil(runiform(K, 1) * min((9, K)))))
    }
    return(phdid_relabel(ceil(runiform(K, 1) * min((3, K)))))
}

// -------------------------------------------------------------------------
// The sampler itself, separated from the reporting so that the
// regularisation paths in phdid_post.mata can run it repeatedly without
// touching Stata's return list.
// -------------------------------------------------------------------------

struct phdid_draws {
    real matrix tau
    real matrix z
    real colvector m
    real colvector chain
    real scalar n
}

// Data-scaled base measure, matching the paper's applications: mu0 is the
// median of the first-stage effects and sigma0^2 is (10 sd(tau))^2.  A
// missing value in either argument asks for the default.
real rowvector phdid_basemeasure(struct phdid_data scalar d, real scalar mu0in,
                                 real scalar s0sqin)
{
    real scalar mu0, s0sq, sprd

    mu0 = mu0in
    if (mu0 >= .) {
        mu0 = mm_median_phdid(d.tau)
    }
    s0sq = s0sqin
    if (s0sq >= .) {
        sprd = sqrt(quadcross(d.tau :- mean(d.tau), d.tau :- mean(d.tau)) / (d.K - 1))
        if (sprd <= 0 | sprd >= .) {
            sprd = max(sqrt(diagonal(d.Sigma)))
        }
        s0sq = (10 * sprd)^2
    }
    return((mu0, s0sq))
}

struct phdid_draws scalar phdid_gibbs(struct phdid_data scalar d,
                                      real scalar alpha, real scalar mu0,
                                      real scalar s0sq, real scalar iters,
                                      real scalar burn, real scalar thin,
                                      real scalar chains, real scalar exact,
                                      string scalar init)
{
    struct phdid_draws scalar D
    real colvector Otau, prec, z, phi, taui
    real scalar K, nkeep, kept, ch, it

    K    = d.K
    Otau = d.Omega * d.tau
    prec = diagonal(d.Omega)

    nkeep = 0
    for (it = burn + 1; it <= iters; it++) {
        if (mod(it - burn - 1, thin) == 0) {
            nkeep = nkeep + 1
        }
    }
    nkeep = nkeep * chains

    D.tau   = J(nkeep, K, .)
    D.z     = J(nkeep, K, .)
    D.m     = J(nkeep, 1, .)
    D.chain = J(nkeep, 1, .)

    kept = 0
    for (ch = 1; ch <= chains; ch++) {
        z = phdid_init(init, K, ch)
        for (it = 1; it <= iters; it++) {
            z = phdid_sweep(z, d.Omega, Otau, d.tau, prec, alpha, mu0, s0sq, exact)
            phi  = phdid_draw_phi(z, d.Omega, Otau, mu0, s0sq)
            taui = phdid_indicator(z) * phi
            if (it > burn & mod(it - burn - 1, thin) == 0) {
                kept = kept + 1
                D.tau[kept, .] = taui'
                D.z[kept, .]   = z'
                D.m[kept]      = max(z)
                D.chain[kept]  = ch
            }
        }
    }
    D.n = kept
    return(D)
}

void phdid_bayes_run(real scalar alpha, real scalar mu0in, real scalar s0sqin,
                     real scalar iters, real scalar burn, real scalar thin,
                     real scalar chains, string scalar marginal,
                     string scalar init, real scalar level)
{
    external real matrix PHDID_TAUDRAWS
    external real matrix PHDID_ZDRAWS
    external real colvector PHDID_MDRAWS
    external real colvector PHDID_CHAIN

    struct phdid_data scalar d
    struct phdid_draws scalar D
    real colvector mdraws, chainid, lo, hi, sdv, qq
    real matrix taudraws, zdraws, cc, res
    real rowvector base
    real scalar K, mu0, s0sq, exact, kept, k, j, a, i

    d = phdid_getdata()
    K = d.K

    base  = phdid_basemeasure(d, mu0in, s0sqin)
    mu0   = base[1]
    s0sq  = base[2]
    exact = (marginal == "exact")

    D = phdid_gibbs(d, alpha, mu0, s0sq, iters, burn, thin, chains, exact, init)
    taudraws = D.tau
    zdraws   = D.z
    mdraws   = D.m
    chainid  = D.chain
    kept     = D.n

    // co-clustering matrix, eq. (39)
    cc = J(K, K, 0)
    for (i = 1; i <= kept; i++) {
        for (j = 1; j <= K; j++) {
            for (k = 1; k <= K; k++) {
                if (zdraws[i, j] == zdraws[i, k]) {
                    cc[j, k] = cc[j, k] + 1
                }
            }
        }
    }
    cc = cc / kept

    a = (1 - level / 100) / 2
    lo  = J(K, 1, .)
    hi  = J(K, 1, .)
    sdv = J(K, 1, .)
    for (k = 1; k <= K; k++) {
        qq = sort(taudraws[., k], 1)
        lo[k] = qq[max((1, ceil(a * kept)))]
        hi[k] = qq[min((kept, ceil((1 - a) * kept)))]
        sdv[k] = sqrt(quadcross(taudraws[., k] :- mean(taudraws[., k]),
                                taudraws[., k] :- mean(taudraws[., k])) / (kept - 1))
    }

    res = J(K, 5, .)
    for (k = 1; k <= K; k++) {
        res[k, 1] = d.tau[k]
        res[k, 2] = mean(taudraws[., k])
        res[k, 3] = sdv[k]
        res[k, 4] = lo[k]
        res[k, 5] = hi[k]
    }

    PHDID_TAUDRAWS = taudraws
    PHDID_ZDRAWS   = zdraws
    PHDID_MDRAWS   = mdraws
    PHDID_CHAIN    = chainid

    st_matrix("r(cells)", res)
    st_matrixcolstripe("r(cells)", (J(5,1,""),
        ("flexible" \ "posterior" \ "sd" \ "lower" \ "upper")))
    st_matrixrowstripe("r(cells)", (J(K,1,""), d.lab))
    st_matrix("r(coclust)", cc)
    st_matrixcolstripe("r(coclust)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(coclust)", (J(K,1,""), d.lab))
    // posterior mean and covariance of the cell effects, posted so that
    // lincom, test and nlcom work on the partition-averaged posterior
    st_matrix("r(b)", mean(taudraws))
    st_matrixcolstripe("r(b)", (J(K,1,""), d.lab))
    st_matrix("r(V)", variance(taudraws))
    st_matrixcolstripe("r(V)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(V)", (J(K,1,""), d.lab))
    st_numscalar("r(m_mean)", mean(mdraws))
    st_numscalar("r(ndraws)", kept)
    st_numscalar("r(K)", K)
    st_numscalar("r(alpha)", alpha)
    st_numscalar("r(mu0)", mu0)
    st_numscalar("r(s0sq)", s0sq)
    st_numscalar("r(prior_m)", alpha * ln(1 + K / alpha))
}

end

mata:

// =====================================================================
// EXACT PARTITION POSTERIOR BY ENUMERATION  (paper Appendix E)
//
// When K is small the posterior can be computed in closed form rather
// than sampled: evaluate the collapsed marginal likelihood and the
// Chinese Restaurant Process prior at every one of the B_K partitions
// and normalise.  This is the deterministic oracle the sampler is
// validated against, and the reason a different RNG does not matter.
// =====================================================================

real scalar phdid_bell(real scalar K)
{
    real colvector row, nw
    real scalar i, j

    row = 1
    for (i = 1; i <= K; i++) {
        nw = J(i + 1, 1, 0)
        nw[1] = row[rows(row)]
        for (j = 1; j <= i; j++) {
            nw[j + 1] = nw[j] + row[j]
        }
        row = nw
    }
    return(row[1])
}

// All set partitions as restricted growth strings: z[1]=1 and
// z[j] <= max(z[1..j-1]) + 1, which enumerates each exactly once.
real matrix phdid_allparts(real scalar K)
{
    real matrix out
    real colvector z, mx
    real scalar n, row, i, done

    n = phdid_bell(K)
    out = J(n, K, 0)
    z  = J(K, 1, 1)
    mx = J(K, 1, 1)
    row = 0

    while (1) {
        row = row + 1
        out[row, .] = z'
        if (row >= n) {
            break
        }
        // advance to the next restricted growth string
        i = K
        done = 0
        while (i >= 2 & done == 0) {
            if (z[i] <= mx[i - 1]) {
                z[i] = z[i] + 1
                done = 1
            }
            else {
                z[i] = 1
                i = i - 1
            }
        }
        for (i = 2; i <= K; i++) {
            if (z[i] > mx[i - 1]) {
                mx[i] = z[i]
            }
            else {
                mx[i] = mx[i - 1]
            }
        }
    }
    return(out)
}

// log CRP prior, eq. (25):
//   Pr(P) = alpha^m prod_l (|C_l|-1)! / alpha^(K)
real scalar phdid_logcrp(real colvector z, real scalar alpha)
{
    real colvector sz
    real scalar m, K, s, i

    sz = phdid_sizes(z)
    m  = rows(sz)
    K  = rows(z)
    s  = m * ln(alpha)
    for (i = 1; i <= m; i++) {
        s = s + lngamma(sz[i])
    }
    for (i = 1; i <= K; i++) {
        s = s - ln(alpha + i - 1)
    }
    return(s)
}

void phdid_enum_run(real scalar alpha, real scalar mu0in, real scalar s0sqin,
                    real scalar maxcells, real scalar level)
{
    struct phdid_data scalar d
    real matrix Z, cc, Sstar, R, S
    real colvector Otau, lp, prob, mvec, aggmu, aggsd, taum, u, mustar, Rw
    real scalar K, n, i, j, mu0, s0sq, sprd, mx, tot, lo, hi, amean, k

    d = phdid_getdata()
    K = d.K
    if (K > maxcells) {
        errprintf("exact enumeration needs all %f partitions of %f cells\n",
                  phdid_bell(K), K)
        errprintf("{p 4 4 2}That is beyond the maxcells() limit. Use ")
        errprintf("{bf:phdid bayes} and check convergence instead.{p_end}\n")
        exit(error(198))
    }

    Otau = d.Omega * d.tau
    mu0 = mu0in
    if (mu0 >= .) {
        mu0 = mm_median_phdid(d.tau)
    }
    s0sq = s0sqin
    if (s0sq >= .) {
        sprd = sqrt(quadcross(d.tau :- mean(d.tau), d.tau :- mean(d.tau)) / (K - 1))
        s0sq = (10 * sprd)^2
    }

    Z = phdid_allparts(K)
    n = rows(Z)
    lp    = J(n, 1, .)
    mvec  = J(n, 1, .)
    aggmu = J(n, 1, .)
    aggsd = J(n, 1, .)
    taum  = J(n, K, .)

    for (i = 1; i <= n; i++) {
        lp[i]   = phdid_logml(Z[i, .]', d.Omega, Otau, mu0, s0sq) +
                  phdid_logcrp(Z[i, .]', alpha)
        mvec[i] = max(Z[i, .])

        R = phdid_indicator(Z[i, .]')
        S = quadcross(R, d.Omega * R)
        u = quadcross(R, Otau)
        Sstar = phdid_spdinv(S + I(cols(R)) / s0sq, "group-effect precision")
        mustar = Sstar * (u :+ mu0 / s0sq)
        taum[i, .] = (R * mustar)'
        aggmu[i] = quadcross(d.w, taum[i, .]')
        Rw = quadcross(R, d.w)
        aggsd[i] = sqrt(max((0, quadcross(Rw, Sstar * Rw))))
    }

    mx = max(lp)
    prob = exp(lp :- mx)
    tot = colsum(prob)
    prob = prob / tot

    cc = J(K, K, 0)
    for (i = 1; i <= n; i++) {
        if (prob[i] < 1e-12) {
            continue
        }
        for (k = 1; k <= K; k++) {
            for (j = 1; j <= K; j++) {
                if (Z[i, k] == Z[i, j]) {
                    cc[k, j] = cc[k, j] + prob[i]
                }
            }
        }
    }
    st_matrix("r(coclust)", cc)
    st_matrixcolstripe("r(coclust)", (J(K,1,""), d.lab))
    st_matrixrowstripe("r(coclust)", (J(K,1,""), d.lab))

    st_numscalar("r(n_partitions)", n)
    st_numscalar("r(m_mean)", quadcross(prob, mvec))
    amean = quadcross(prob, aggmu)
    st_numscalar("r(att)", amean)

    // the aggregate posterior is a finite mixture of normals, one per
    // partition; its quantiles come from inverting the mixture CDF
    lo = phdid_mixq(prob, aggmu, aggsd, (1 - level / 100) / 2)
    hi = phdid_mixq(prob, aggmu, aggsd, 1 - (1 - level / 100) / 2)
    st_numscalar("r(att_lo)", lo)
    st_numscalar("r(att_hi)", hi)
    st_matrix("r(prob)", prob')
    st_matrix("r(tau)", (prob' * taum))
    st_numscalar("r(K)", K)
    st_numscalar("r(alpha)", alpha)
}

real scalar phdid_mixq(real colvector prob, real colvector mu,
                       real colvector sd, real scalar q)
{
    real scalar lo, hi, mid, i, cdf, sdm

    sdm = max((1e-12, max(sd)))
    lo = min(mu) - 10 * sdm
    hi = max(mu) + 10 * sdm
    for (i = 1; i <= 200; i++) {
        mid = (lo + hi) / 2
        cdf = colsum(prob :* normal((mid :- mu) :/ (sd :+ 1e-300)))
        if (cdf < q) {
            lo = mid
        }
        else {
            hi = mid
        }
    }
    return((lo + hi) / 2)
}

end
