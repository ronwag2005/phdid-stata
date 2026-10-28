*! phdid_core.mata -- numerical core for the phdid package
*! Partial homogeneity in staggered difference-in-differences
*! Rohan Wagle and Parush Arora
version 17.0

mata: mata set matastrict on

mata:

// =====================================================================
// SESSION PROBE
//
// phdid_mata.ado calls this to find out whether the library is already
// compiled in the session.  The number is an interface generation, to
// be bumped only if a future release changes a function signature that
// the ado-files call.
// =====================================================================

real scalar phdid_abi()
{
    return(1)
}

// =====================================================================
// PARTITION HANDLING
//
// A partition of the K cohort-time cells is an integer column vector z
// of length K taking values in 1..m, canonicalised so that labels
// appear in order of first occurrence.  The group-indicator matrix R is
// K x m with R[k,p] = 1 when cell k belongs to group p; partial
// homogeneity is the restriction tau = R*phi.
// =====================================================================

real colvector phdid_relabel(real colvector z)
{
    real colvector out, seen
    real scalar i, nxt, K

    K   = rows(z)
    out = J(K, 1, 0)
    seen = J(max(z), 1, 0)
    nxt = 0
    for (i = 1; i <= K; i++) {
        if (seen[z[i]] == 0) {
            nxt = nxt + 1
            seen[z[i]] = nxt
        }
        out[i] = seen[z[i]]
    }
    return(out)
}

real matrix phdid_indicator(real colvector z)
{
    real colvector zz
    real matrix R
    real scalar i

    zz = phdid_relabel(z)
    R  = J(rows(zz), max(zz), 0)
    for (i = 1; i <= rows(zz); i++) {
        R[i, zz[i]] = 1
    }
    return(R)
}

real colvector phdid_sizes(real colvector z)
{
    real colvector zz, sz
    real scalar i

    zz = phdid_relabel(z)
    sz = J(max(zz), 1, 0)
    for (i = 1; i <= rows(zz); i++) {
        sz[zz[i]] = sz[zz[i]] + 1
    }
    return(sz)
}

// Number of cross-group pairs c(P) = C(K,2) - sum_p C(|C_p|,2),
// the complexity term of the l0 penalty (paper eq. 33).
real scalar phdid_crosspairs(real colvector z)
{
    real colvector sz
    real scalar p, tot, K

    K   = rows(z)
    sz  = phdid_sizes(z)
    tot = 0
    for (p = 1; p <= rows(sz); p++) {
        tot = tot + sz[p] * (sz[p] - 1) / 2
    }
    return(K * (K - 1) / 2 - tot)
}

string scalar phdid_fmtpart(real colvector z, string colvector lab)
{
    real colvector zz
    real scalar p, k, first
    string scalar s

    zz = phdid_relabel(z)
    s  = ""
    for (p = 1; p <= max(zz); p++) {
        if (p > 1) {
            s = s + " "
        }
        s = s + "{"
        first = 1
        for (k = 1; k <= rows(zz); k++) {
            if (zz[k] == p) {
                if (first == 0) {
                    s = s + ","
                }
                s = s + lab[k]
                first = 0
            }
        }
        s = s + "}"
    }
    return(s)
}

// =====================================================================
// RESTRICTED GLS UNDER A PARTITION  (paper eq. 12, 50, 55)
//
//   S    = R' Omega R            m x m grouped precision
//   u    = R' Omega tau          m x 1 grouped score
//   phi  = S^-1 u                group effects
//   tau~ = R phi                 cell effects
//   V    = R S^-1 R'             cell covariance
//   dev  = tau'Omega tau - u'phi GLS deviance
//
// This is the workhorse.  Averaging the flexible estimates within each
// group coincides with it only when the within-transformed dummies are
// orthogonal; in a general panel the group means fail the normal
// equations and are not best linear unbiased (paper Appendix B).
// =====================================================================

struct phdid_fit {
    real colvector tau
    real colvector phi
    real colvector z
    real matrix V
    real matrix Vphi
    real scalar m
    real scalar dev
}

real matrix phdid_spdinv(real matrix A, string scalar what)
{
    real matrix Asym, Ainv

    // cholinv() takes the symmetric positive-definite matrix itself, not
    // its Cholesky factor, and returns missing values when A is not
    // positive definite.
    Asym = (A + A') / 2
    Ainv = cholinv(Asym)
    if (missing(Ainv) > 0) {
        errprintf("the %s is not positive definite\n", what)
        errprintf("{p 4 4 2}This usually means the first-stage covariance ")
        errprintf("is singular, which happens when a cohort-time cell has ")
        errprintf("no clean comparison group and so carries no independent ")
        errprintf("identifying variation. Drop the affected cells, or ")
        errprintf("supply a regularised covariance.{p_end}\n")
        exit(error(506))
    }
    return(Ainv)
}

struct phdid_fit scalar phdid_gls(real colvector tau, real matrix Omega,
                                  real colvector z)
{
    struct phdid_fit scalar f
    real matrix R, S, Sinv, OR
    real colvector u

    f.z  = phdid_relabel(z)
    R    = phdid_indicator(f.z)
    OR   = Omega * R
    S    = quadcross(R, OR)
    u    = quadcross(OR, tau)
    Sinv = phdid_spdinv(S, "grouped precision R'Sigma^-1 R")

    f.phi  = Sinv * u
    f.tau  = R * f.phi
    f.Vphi = Sinv
    f.V    = R * Sinv * R'
    f.m    = cols(R)
    f.dev  = quadcross(tau, Omega * tau) - quadcross(u, f.phi)
    return(f)
}

// Information criterion for a partition.  On the two-stage route the
// first-stage covariance is taken as known, so the scale is fixed and
// the deviance form applies: BIC = dev + m*ln(n).  The sample size n is
// not forced by the theory -- the K first-stage estimates are the
// observations, while the asymptotics justifying Schwarz's rate run in
// the number of units -- so it is exposed rather than buried.  K is the
// default, which reproduces the paper's reported group counts.
real scalar phdid_bic(real scalar dev, real scalar m, real scalar nbic)
{
    return(dev + m * ln(nbic))
}

// =====================================================================
// L0-PENALISED ESTIMATOR  (paper Propositions 2-3, Appendix B)
//
// Step 1  greedy agglomerative search from K singletons, merging at each
//         stage the pair minimising
//             dObj(A,B) = n_A n_B / (n_A + n_B) * (tau_A - tau_B)^2
// Step 2  exact GLS re-estimation under the selected partition.
//
// Step 2 is not a refinement: the merge criterion of step 1 treats the
// design as orthogonal, and only the re-estimated fit is BLUE in a
// general panel.
// =====================================================================

struct phdid_path {
    real matrix parts
    real colvector mergecost
    real colvector pairprod
    real colvector effsize
}

struct phdid_path scalar phdid_agglom(real colvector tau, real matrix Omega,
                                      real scalar sigma2)
{
    struct phdid_path scalar P
    real colvector gtau, gn, gsize, prec, z
    real matrix members
    real scalar K, step, ng, a, b, best, bi, bj, cost, newn, i, j

    K    = rows(tau)
    prec = diagonal(Omega)
    if (sigma2 < .) {
        prec = prec * sigma2
    }

    members = I(K)
    gtau    = tau
    gn      = prec
    gsize   = J(K, 1, 1)

    P.parts     = J(K, K, .)
    P.mergecost = J(K, 1, .)
    P.pairprod  = J(K, 1, .)
    P.effsize   = prec
    P.parts[., K] = (1::K)

    ng = K
    for (step = 1; step <= K - 1; step++) {
        best = .
        bi   = .
        bj   = .
        for (a = 1; a <= ng - 1; a++) {
            for (b = a + 1; b <= ng; b++) {
                cost = (gn[a] * gn[b] / (gn[a] + gn[b])) * (gtau[a] - gtau[b])^2
                if (cost < best) {
                    best = cost
                    bi   = a
                    bj   = b
                }
            }
        }
        P.mergecost[ng - 1] = best
        P.pairprod[ng - 1]  = gsize[bi] * gsize[bj]

        newn      = gn[bi] + gn[bj]
        gtau[bi]  = (gn[bi] * gtau[bi] + gn[bj] * gtau[bj]) / newn
        gn[bi]    = newn
        gsize[bi] = gsize[bi] + gsize[bj]
        members[., bi] = members[., bi] + members[., bj]

        if (bj < ng) {
            for (j = bj; j <= ng - 1; j++) {
                members[., j] = members[., j + 1]
                gtau[j]  = gtau[j + 1]
                gn[j]    = gn[j + 1]
                gsize[j] = gsize[j + 1]
            }
        }
        ng = ng - 1

        z = J(K, 1, 0)
        for (j = 1; j <= ng; j++) {
            for (i = 1; i <= K; i++) {
                if (members[i, j] > 0) {
                    z[i] = j
                }
            }
        }
        P.parts[., ng] = phdid_relabel(z)
    }
    return(P)
}

// Select a point on the path by the Appendix B stopping rule at a fixed
// lambda: merging continues while a pair's cost falls below
// lambda*|A|*|B|.
//
// The effective per-merge threshold cost/(|A||B|) is NOT monotone in the
// merge order, because |A||B| grows as groups absorb one another.  When
// a later merge is cheaper in threshold terms than an earlier one,
// raising lambda past the earlier merge triggers both at once, and no
// value of lambda selects the partition in between.  Index the path by
// m to visit every group count.
real scalar phdid_select_lambda(struct phdid_path scalar P,
                                real scalar lambda, real scalar K)
{
    real scalar m, i

    m = K
    for (i = K - 1; i >= 1; i--) {
        if (P.mergecost[i] >= .) {
            break
        }
        if (P.mergecost[i] < lambda * P.pairprod[i]) {
            m = i
        }
        else {
            break
        }
    }
    return(m)
}

end
