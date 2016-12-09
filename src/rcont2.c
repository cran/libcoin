
#include "libcoin_internal.h"

#ifdef USE_RCONT2_FROM_R
/* This R _header_ file containts the function defintion for S_rcont2, a wrapper
   for rcont2 in stats/src. We need to compile this for libcoin. Looks like this   
   mechanism was invented for lme4 < 1.0-0 */ 
#include <R_ext/stats_stubs.h>
#else
/* The above export is not there yet, so just copy rcont2 */
void S_rcont2
(
    int *nrow,
    int *ncol,
    int *nrowt,
    int *ncolt,
    int *ntotal,
    double *fact,
    int *jwork,
    int *matrix
) {

    int j, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm, nr_1, nc_1;
    double x, y, dummy, sumprb;
    Rboolean lsm, lsp;

    nr_1 = *nrow - 1;
    nc_1 = *ncol - 1;

    ib = 0; /* -Wall */

    /* Construct random matrix */
    for (j = 0; j < nc_1; ++j)
        jwork[j] = ncolt[j];

    jc = *ntotal;

    for (l = 0; l < nr_1; ++l) { /* -----  matrix[ l, * ] ----- */
        ia = nrowt[l];
        ic = jc;
        jc -= ia;/* = n_tot - sum(nr[0:l]) */

        for (m = 0; m < nc_1; ++m) {
            id = jwork[m];
            ie = ic;
            ic -= id;
            ib = ie - ia;
            ii = ib - id;

            if (ie == 0) { /* Row [l,] is full, fill rest with zero entries */
                for (j = m; j < nc_1; ++j)
                    matrix[l + j * *nrow] = 0;
                ia = 0;
                break;
            }

            /* Generate pseudo-random number */
            dummy = unif_rand();

            do {/* Outer Loop */

                /* Compute conditional expected value of MATRIX(L, M) */

                nlm = (int)(ia * (id / (double) ie) + 0.5);
                x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]
                        - fact[ie] - fact[nlm]
                        - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);
                if (x >= dummy)
                    break;
                if (x == 0.)/* MM: I haven't seen this anymore */
                    error("rcont2 [%d,%d]: exp underflow to 0; algorithm failure", l, m);

                sumprb = x;
                y = x;
                nll = nlm;

                do {
                    /* Increment entry in row L, column M */
                    j = (int)((id - nlm) * (double)(ia - nlm));
                    lsp = (j == 0);
                    if (!lsp) {
                        ++nlm;
                        x = x * j / ((double) nlm * (ii + nlm));
                        sumprb += x;
                        if (sumprb >= dummy)
                            goto L160;
                    }

                    do {
                        R_CheckUserInterrupt();

                        /* Decrement entry in row L, column M */
                        j = (int)(nll * (double)(ii + nll));
                        lsm = (j == 0);
                        if (!lsm) {
                            --nll;
                            y = y * j / ((double) (id - nll) * (ia - nll));
                            sumprb += y;
                            if (sumprb >= dummy) {
                                nlm = nll;
                                goto L160;
                            }
                            /* else */
                            if (!lsp)
                                break;/* to while (!lsp) */
                        }
                    } while (!lsm);

                } while (!lsp);

                dummy = sumprb * unif_rand();

            } while (1);

L160:
            matrix[l + m * *nrow] = nlm;
            ia -= nlm;
            jwork[m] -= nlm;
        }
        matrix[l + nc_1 * *nrow] = ia;/* last column in row l */
    }

    /* Compute entries in last row of MATRIX */
    for (m = 0; m < nc_1; ++m)
        matrix[nr_1 + m * *nrow] = jwork[m];

    matrix[nr_1 + nc_1 * *nrow] = ib - matrix[nr_1 + (nc_1-1) * *nrow];

    return;
}
#endif
