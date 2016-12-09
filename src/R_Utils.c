
#include "libcoin_internal.h"
#include "Utils.h"

/**
    R-interface to C_kronecker\n
    *\param A matrix
    *\param B matrix
*/

SEXP R_kronecker
(
    SEXP A,
    SEXP B
) {

    int m, n, r, s;
    SEXP ans;

    if (!isReal(A) || !isReal(B))
        error("R_kronecker: A and / or B are not of type REALSXP");

    m = NROW(A);
    n = NCOL(A);
    r = NROW(B);
    s = NCOL(B);

    PROTECT(ans = allocVector(REALSXP, m * n * r * s));
    C_kronecker(REAL(A), m, n, REAL(B), r, s, 1, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
