
#include "libcoin_internal.h"
#include "LinearStatistic.h"
#include "TestStatistics.h"
#include "Utils.h"
#include "Sums.h"
#include "Tables.h"
#include "Distributions.h"
#include "MemoryAccess.h"
#include <R_ext/stats_stubs.h> /* for S_rcont2 */

void RC_ExpectationCovarianceStatistic
(
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    SEXP ans
) {

    int N, P, Q, Lb, *sumweights, *table, *subset_tmp, tmp;
    double *ExpInf, *work;

    /* note: x being an integer (Xfactor) with some 0 elements is not
             handled correctly (as sumweights doesnt't take this information
             into account; use subset to exclude these missings (as done
             in libcoin::LinStatExpCov) */

    P = C_get_P(ans);
    Q = C_get_Q(ans);

    N = NROW(x);
    Lb = 1;
    if (LENGTH(block) > 0)
        Lb = NLEVELS(block);

    ExpInf = C_get_ExpectationInfluence(ans);
    work = C_get_Work(ans);
    table = C_get_TableBlock(ans);
    sumweights = C_get_Sumweights(ans);

    if (Lb == 1) {
        table[0] = 0;
        table[1] = -1; /* means: no subset given */
        if (LENGTH(subset) > 0)
            table[1] = LENGTH(subset);
        if (LENGTH(weights) == 0) {
            sumweights[0] = -1; /* means: no weights given */
        } else {
            if (LENGTH(subset) == 0) {
                sumweights[0] = C_sum_weights(INTEGER(weights), LENGTH(weights));
            } else {
                sumweights[0] = C_sum_weights_subset(INTEGER(weights), LENGTH(weights),
                                             INTEGER(subset), LENGTH(subset));
            }
        }
        subset_tmp = INTEGER(subset);
    } else {
        if (LENGTH(subset) == 0) {
            C_1dtable_(INTEGER(block), Lb + 1, N, table);
            subset_tmp = Calloc(N, int);
            C_setup_subset(N, subset_tmp);
            C_order_wrt_block(subset_tmp, N, INTEGER(block), table, Lb + 1);
        } else {
            C_1dtable_subset(INTEGER(block), Lb + 1, INTEGER(subset),
                             LENGTH(subset), table);
            subset_tmp = Calloc(LENGTH(subset), int);
            Memcpy(subset_tmp, INTEGER(subset), LENGTH(subset));
            C_order_wrt_block(subset_tmp, LENGTH(subset), INTEGER(block),
                              table, Lb + 1);
        }

        if (LENGTH(weights) == 0) {
            for (int b = 0; b < Lb; b++) sumweights[b] = -1; /* means: no weights given */
        } else {
            tmp = 0;
            for (int b = 0; b < Lb; b++) {
                sumweights[b] = C_sum_weights_subset(INTEGER(weights), LENGTH(weights),
                                             subset_tmp + tmp, table[b + 1]);
                tmp = tmp + table[b + 1];
            }
        }
    }

    RC_LinearStatistic(x, N, P, REAL(y), Q, INTEGER(weights),
                      sumweights, subset_tmp, table + 1, Lb, C_get_LinearStatistic(ans));

    C_ExpectationCoVarianceInfluence(REAL(y), N, Q, INTEGER(weights),
                                     sumweights, subset_tmp, table + 1, Lb, 1, ExpInf,
                                     C_get_VarianceInfluence(ans), C_get_CovarianceInfluence(ans));

    if (C_get_varonly(ans)) {
        RC_ExpectationVarianceLinearStatistic(x, N, P, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Lb, C_get_ExpectationX(ans), ExpInf, C_get_VarianceInfluence(ans),
            work, C_get_Expectation(ans), C_get_Variance(ans));
    } else {
        RC_ExpectationCovarianceLinearStatistic(x, N, P, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Lb, C_get_ExpectationX(ans), ExpInf,
            C_get_CovarianceInfluence(ans), work,
            C_get_Expectation(ans), C_get_Covariance(ans));
    }
    if (Lb > 1) Free(subset_tmp);
}


SEXP R_ExpectationCovarianceStatistic
(
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP varonly,
    const SEXP tol
) {

    SEXP ans, P, Q, Lb, Xfactor;
    R_xlen_t n;

    /* check for long vectors (not supported at the moment) */
    n = xlength(y);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");
    n = xlength(subset);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");

    PROTECT(P = ScalarInteger(0));
    PROTECT(Q = ScalarInteger(0));
    PROTECT(Lb = ScalarInteger(0));
    PROTECT(Xfactor = ScalarInteger(0));

    if (isInteger(x)) {
        INTEGER(P)[0] = NLEVELS(x);
        INTEGER(Xfactor)[0] = 1;
    } else {
        INTEGER(P)[0] = NCOL(x);
        INTEGER(Xfactor)[0] = 0;
    }
    INTEGER(Q)[0] = NCOL(y);

    INTEGER(Lb)[0] = 1;
    if (LENGTH(block) > 0)
        INTEGER(Lb)[0] = NLEVELS(block);

    PROTECT(ans = R_init_LECV_1d(P, Q, varonly, Lb, Xfactor, tol));

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, ans);

    UNPROTECT(5);
    return(ans);
}

SEXP R_PermutedLinearStatistic
(
    const SEXP LEV,
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP B,
    const SEXP standardise
) {

    SEXP ans;
    double *linstat;
    int *orig, *perm, *tmp;
    int P, Q, PQ, Lb, N, *table;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
    Lb = 1;
    if (LENGTH(block) > 0)
        Lb = NLEVELS(block);

    PROTECT(ans = allocMatrix(REALSXP, PQ, INTEGER(B)[0]));

    GetRNGstate();

    if (Lb == 1) {
        if (LENGTH(weights) == 0) {
            if (LENGTH(subset) == 0) {
                N = NROW(x);
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_setup_subset(N, orig);
            } else {
                N = LENGTH(subset);
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                Memcpy(orig, INTEGER(subset), N);
            }
        } else {
            if (LENGTH(subset) == 0) {
                N = C_sum_weights(INTEGER(weights), LENGTH(weights));
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_setup_subset_weights(LENGTH(weights), INTEGER(weights), orig);
            } else {
                N = C_sum_weights_subset(INTEGER(weights), LENGTH(weights),
                                 INTEGER(subset), LENGTH(subset));
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_setup_subset_weights_subset(LENGTH(subset), INTEGER(weights),
                                              INTEGER(subset), orig);
            }
        }
        tmp = Calloc(N, int);
        for (int i = 0; i < INTEGER(B)[0]; i++) {
            if (i % 256 == 0) R_CheckUserInterrupt();
            linstat = REAL(ans) + PQ * i;
            for (int p = 0; p < PQ; p++)
                linstat[p] = 0;
            C_doPermute(orig, N, tmp, perm);
            RC_PermutedLinearStatistic(x, NROW(x), P, REAL(y), Q, perm, orig, N, linstat);
        }
    } else {
        table = Calloc(Lb + 1, int);
        if (LENGTH(weights) == 0) {
            if (LENGTH(subset) == 0) {
                N = LENGTH(block);
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_1dtable_(INTEGER(block), Lb + 1, LENGTH(block), table);
                C_setup_subset(N, orig);
                C_order_wrt_block(orig, N, INTEGER(block), table, Lb + 1);
            } else {
                N = LENGTH(subset);
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_1dtable_subset(INTEGER(block), Lb + 1, INTEGER(subset), N,
                                 table);
                Memcpy(orig, INTEGER(subset), N);
                C_order_wrt_block(orig, N, INTEGER(block), table, Lb + 1);
            }
        } else {
            if (LENGTH(subset) == 0) {
                N = C_sum_weights(INTEGER(weights), LENGTH(weights));
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_1dtable_weights(INTEGER(block), Lb + 1, INTEGER(weights),
                                  LENGTH(weights), table);
                C_setup_subset_weights(LENGTH(weights), INTEGER(weights),
                                       orig);
                C_order_wrt_block(orig, N, INTEGER(block), table, Lb + 1);
            } else {
                N = C_sum_weights_subset(INTEGER(weights), LENGTH(weights),
                                 INTEGER(subset), LENGTH(subset));
                orig = Calloc(N, int);
                perm = Calloc(N, int);
                C_1dtable_weights_subset(INTEGER(block), Lb + 1,
                                         INTEGER(weights), INTEGER(subset),
                                         LENGTH(subset), table);
                C_setup_subset_weights_subset(LENGTH(subset), INTEGER(weights),
                                              INTEGER(subset), orig);
                C_order_wrt_block(orig, N, INTEGER(block), table, Lb + 1);
            }
        }
        tmp = Calloc(N, int);
        for (int i = 0; i < INTEGER(B)[0]; i++) {
            if (i % 256 == 0) R_CheckUserInterrupt();
            linstat = REAL(ans) + PQ * i;
            for (int p = 0; p < PQ; p++)
                linstat[p] = 0;
            C_doPermuteBlock(orig, N, table, Lb + 1, tmp, perm);
            RC_PermutedLinearStatistic(x, NROW(x), P, REAL(y), Q,
                                      perm, orig, N, linstat);
        }
        Free(table);
    }
    PutRNGstate();
    Free(tmp); Free(perm); Free(orig);

    if (INTEGER(standardise)[0]) {
        for (int i = 0; i < INTEGER(B)[0]; i++) {
            if (C_get_varonly(LEV)) {
                C_standardise(PQ, REAL(ans) + PQ * i, C_get_Expectation(LEV),
                              C_get_Variance(LEV), 1, C_get_tol(LEV));
            } else {
                C_standardise(PQ, REAL(ans) + PQ * i, C_get_Expectation(LEV),
                              C_get_Covariance(LEV), 0, C_get_tol(LEV));
            }
        }
    }

    UNPROTECT(1);
    return(ans);
}

void RC_ExpectationCovarianceStatistic_2d
(
    const SEXP x,
    const SEXP ix,
    const SEXP y,
    const SEXP iy,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    SEXP ans
) {

    int P, Q, Lxp1, Lyp1, Lb, *btab, *csum, *rsum, *table, *table2d, sw;
    double *ExpInf, *ExpX, *CovX, *work;


    P = C_get_P(ans);
    Q = C_get_Q(ans);

    ExpX = C_get_ExpectationX(ans);
    table = C_get_Table(ans);

    Lxp1 = C_get_dimTable(ans)[0];
    Lyp1 = C_get_dimTable(ans)[1];
    Lb = C_get_dimTable(ans)[2];

    table2d = Calloc(Lxp1 * Lyp1, int);
    csum = Calloc(Lyp1, int);
    rsum = Calloc(Lxp1, int);

    for (int i = 0; i < Lxp1 * Lyp1; i++)
        table2d[i] = 0;
    for (int b = 0; b < Lb; b++) {
        for (int i = 0; i < Lxp1; i++) {
            for (int j = 0; j < Lyp1; j++)
                table2d[j * Lxp1 + i] += table[b * Lxp1 * Lyp1 + j * Lxp1 + i];
        }
    }

    ExpInf = C_get_ExpectationInfluence(ans);
    work = C_get_Work(ans);

    if (C_get_varonly(ans)) {
        CovX = work + P;
    } else {
        CovX = work + P * (P + 1) / 2;
    }

    if (LENGTH(x) == 0) {
        RC_LinearStatistic_2d(ix, LENGTH(ix), P, REAL(y), NROW(y), Q,
                              table2d, C_get_LinearStatistic(ans));
    } else {
        RC_LinearStatistic_2d(x, NROW(x), P, REAL(y), NROW(y), Q,
                              table2d, C_get_LinearStatistic(ans));
    }

    for (int b = 0; b < Lb; b++) {
        btab = table + Lxp1 * Lyp1 * b;
        C_colSums_i(btab, Lxp1, Lyp1, csum);
        csum[0] = 0; /* NA */
        C_rowSums_i(btab, Lxp1, Lyp1, rsum);
        rsum[0] = 0; /* NA */
        sw = 0;
        for (int i = 1; i < Lxp1; i++) sw += rsum[i];
        C_ExpectationInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf);
        if (LENGTH(x) == 0) {
            for (int p = 0; p < P; p++)
                ExpX[p] = (double) rsum[p + 1];
        } else {
            C_ExpectationX_weights(REAL(x), NROW(x), P, rsum, ExpX);
        }
        C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, C_get_Expectation(ans));
        /* compute both for the time being */
        C_VarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf,
                                    C_get_VarianceInfluence(ans));
        C_CovarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf,
                                      C_get_CovarianceInfluence(ans));
        if (C_get_varonly(ans)) {
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P; p++) CovX[p] = ExpX[p];
            } else {
                C_VarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            }
            C_VarianceLinearStatistic(P, Q, C_get_VarianceInfluence(ans),
                                      ExpX, CovX, sw, work, b,
                                      C_get_Variance(ans));
        } else {
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P * (P + 1) / 2; p++) CovX[p] = 0.0;
                for (int p = 0; p < P; p++) CovX[S(p, p, P)] = ExpX[p];
            } else {
                C_CovarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            }
            C_CovarianceLinearStatistic(P, Q, C_get_CovarianceInfluence(ans),
                                        ExpX, CovX, sw, work, b,
                                        C_get_Covariance(ans));
        }
    }
    Free(table2d); Free(csum); Free(rsum);
}

SEXP R_ExpectationCovarianceStatistic_2d
(
    const SEXP x,
    const SEXP ix,
    const SEXP y,
    const SEXP iy,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP varonly,
    const SEXP tol
) {

    SEXP ans, P, Q, Lx, Ly, Lb, Xfactor;
    R_xlen_t n;

    /* check for long vectors (not supported at the moment) */
    n = xlength(ix);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");
    n = xlength(subset);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");

    PROTECT(P = ScalarInteger(0));
    PROTECT(Q = ScalarInteger(0));
    PROTECT(Lx = ScalarInteger(0));
    PROTECT(Ly = ScalarInteger(0));
    PROTECT(Lb = ScalarInteger(0));
    PROTECT(Xfactor = ScalarInteger(0));

    if (LENGTH(x) == 0) {
        INTEGER(P)[0] = NLEVELS(ix);
        INTEGER(Xfactor)[0] = 1;
    } else {
        INTEGER(P)[0] = NCOL(x);
        INTEGER(Xfactor)[0] = 0;
    }
    INTEGER(Q)[0] = NCOL(y);

    INTEGER(Lb)[0] = 1;
    if (LENGTH(block) > 0)
        INTEGER(Lb)[0] = NLEVELS(block);

    INTEGER(Lx)[0] = NLEVELS(ix);
    INTEGER(Ly)[0] = NLEVELS(iy);

    PROTECT(ans = R_init_LECV_2d(P, Q, varonly, Lx, Ly, Lb, Xfactor, tol));

    RC_2dtable(ix, iy, weights, subset, block, C_get_Table(ans));
    RC_ExpectationCovarianceStatistic_2d(x, ix, y, iy, weights,
                                         subset, block, ans);

    UNPROTECT(7);
    return(ans);
}

SEXP R_PermutedLinearStatistic_2d
(
    const SEXP LEV,
    const SEXP x,
    const SEXP ix,
    const SEXP y,
    const SEXP iy,
    const SEXP block,
    const SEXP B,
    const SEXP standardise
) {

    SEXP ans;
    int P, Q, PQ, Lb, Lx, Ly, *csum, *rsum, *ntotal, *table, *jwork, *rtable, *rtable2, maxn = 0, Lxp1, Lyp1;
    double *fact, *linstat, *blinstat, *dans;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
    Lxp1 = C_get_dimTable(LEV)[0];
    Lyp1 = C_get_dimTable(LEV)[1];
    Lx = Lxp1 - 1;
    Ly = Lyp1 - 1;
    Lb = C_get_dimTable(LEV)[2];

    table = C_get_Table(LEV);

    PROTECT(ans = allocMatrix(REALSXP, PQ, INTEGER(B)[0]));

    csum = Calloc(Lyp1 * Lb, int);
    rsum = Calloc(Lxp1 * Lb, int);
    ntotal = Calloc(Lb, int);
    rtable = Calloc(Lxp1 * Lyp1, int);
    rtable2 = Calloc(NLEVELS(ix) * NLEVELS(iy) , int);
    linstat = Calloc(PQ, double);
    jwork = Calloc(Lyp1, int);

    for (int b = 0; b < Lb; b++) {
        C_colSums_i(table + Lxp1 * Lyp1 * b,
                    NLEVELS(ix) + 1, NLEVELS(iy) + 1, csum + Lyp1 * b);
        csum[Lyp1 * b] = 0; /* NA */
        C_rowSums_i(table + Lxp1 * Lyp1 * b,
                    Lxp1, Lyp1, rsum + Lxp1 * b);
        rsum[Lxp1 * b] = 0; /* NA */
        ntotal[b] = 0;
        for (int i = 1; i < Lxp1; i++)
            ntotal[b] += rsum[Lxp1 * b + i];
        if (ntotal[b] > maxn) maxn = ntotal[b];
    }

    fact = Calloc(maxn + 1, double);
    /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
    fact[0] = fact[1] = 0.;
    for(int j = 2; j <= maxn; j++)
        fact[j] = fact[j - 1] + log(j);

    GetRNGstate();

    dans = REAL(ans);
    for (int i = 0; i < INTEGER(B)[0]; i++) {

        blinstat = dans + PQ * i;
        for (int p = 0; p < PQ; p++) {
            blinstat[p] = 0.0;
            linstat[p] = 0.0;
        }
        for (int p = 0; p < Lxp1 * Lyp1; p++)
            rtable[p] = 0;

        for (int b = 0; b < Lb; b++) {
            S_rcont2(&Lx, &Ly, rsum + Lxp1 * b + 1,
                     csum + Lyp1 *b + 1, ntotal + b, fact, jwork, rtable2);

        for (int j1 = 1; j1 <= NLEVELS(ix); j1++) {
            for (int j2 = 1; j2 <= NLEVELS(iy); j2++)
                rtable[j2 * Lxp1 + j1] = rtable2[(j2 - 1) * NLEVELS(ix) + (j1 - 1)];
        }
        RC_LinearStatistic_2d(x, Lxp1, P, REAL(y), NROW(y), Q, rtable, linstat);

        for (int p = 0; p < PQ; p++)
            blinstat[p] += linstat[p];
        }
    }

    PutRNGstate();

    if (INTEGER(standardise)[0]) {
        for (int i = 0; i < INTEGER(B)[0]; i++) {
            if (C_get_varonly(LEV)) {
                C_standardise(PQ, REAL(ans) + PQ * i, C_get_Expectation(LEV),
                              C_get_Variance(LEV), 1, C_get_tol(LEV));
            } else {
                C_standardise(PQ, REAL(ans) + PQ * i, C_get_Expectation(LEV),
                              C_get_Covariance(LEV), 0, C_get_tol(LEV));
            }
        }
    }

    Free(csum); Free(rsum); Free(ntotal); Free(rtable); Free(rtable2); Free(linstat);
    Free(jwork); Free(fact);
    UNPROTECT(1);
    return(ans);
}
