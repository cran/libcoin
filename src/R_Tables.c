
#include "libcoin_internal.h"
#include "Tables.h"
#include "Utils.h"

SEXP R_tables
(
    SEXP ix,
    SEXP iy,
    SEXP weights,
    SEXP subset,
    SEXP block
) {

    SEXP ans, tabdim;
    int Lx, Ly, Lb;
    R_xlen_t n;

    /* check for long vectors (not supported at the moment) */
    n = xlength(ix);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");
    n = xlength(subset);
    if (n > INT_MAX)
        error("libcoin does not support long vectors");

    if (LENGTH(iy) > 0) {
        if (LENGTH(ix) != LENGTH(iy))
            error("ix and iy have different length");
    }
    if (LENGTH(weights) > 0) {
        if (LENGTH(ix) != LENGTH(weights))
            error("ix and weights have different length");
    }
    if (LENGTH(block) > 0) {
        if (LENGTH(ix) != LENGTH(block))
            error("ix and block have different length");
    }

    Lx = NLEVELS(ix);
    Ly = 0;
    if (LENGTH(iy) > 0) Ly = NLEVELS(iy);
    Lb = 0;
    if (LENGTH(block) > 0) Lb = NLEVELS(block);

    PROTECT(tabdim = allocVector(INTSXP, 3));
    INTEGER(tabdim)[0] = Lx + 1;
    INTEGER(tabdim)[1] = Ly + 1;
    INTEGER(tabdim)[2] = Lb + 1; /* reuse iy as block */

    if (isInteger(weights)) {
        PROTECT(ans = allocVector(INTSXP, (Lx + 1) * (Ly + 1) * (Lb + 1)));
        dimgets(ans, tabdim);
        if (LENGTH(iy) == 0) {
           RC_1dtable(ix, weights, subset, block, INTEGER(ans));
        } else {
           RC_2dtable(ix, iy, weights, subset, block, INTEGER(ans));
        }
    } else {
        PROTECT(ans = allocVector(REALSXP, (Lx + 1) * (Ly + 1) * (Lb + 1)));
        dimgets(ans, tabdim);
        if (LENGTH(iy) == 0) {
           RC_1dtable_dweights(ix, weights, subset, block, REAL(ans));
        } else {
           RC_2dtable_dweights(ix, iy, weights, subset, block, REAL(ans));
        }
    }

    UNPROTECT(2);
    return(ans);
}
