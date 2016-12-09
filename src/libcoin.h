
#include "libcoin_internal.h"

/* R_LinearStatistic.c */
extern SEXP R_ExpectationCovarianceStatistic
(
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP varonly,
    const SEXP tol
);

extern SEXP R_PermutedLinearStatistic
(
    const SEXP LEV,
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP B,
    const SEXP standardise
);

extern SEXP R_ExpectationCovarianceStatistic_2d
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
);

extern SEXP R_PermutedLinearStatistic_2d
(
    const SEXP LEV,
    const SEXP x,
    const SEXP ix,
    const SEXP y,
    const SEXP iy,
    const SEXP block,
    const SEXP B,
    const SEXP standardise
);


/* R_Tables.c */
extern SEXP R_tables
(
    SEXP ix,
    SEXP iy,
    SEXP weights,
    SEXP subset,
    SEXP block
);


/* R_Tests.c */
extern SEXP R_QuadraticTest
(
    SEXP LEV,
    SEXP pvalue,
    SEXP lower,
    SEXP give_log
);

extern SEXP R_MaximumTest
(
    SEXP LEV,
    SEXP alternative,
    SEXP pvalue,
    SEXP lower,
    SEXP give_log,
    SEXP maxpts,
    SEXP releps,
    SEXP abseps
);

extern SEXP R_MaximallySelectedTest
(
    SEXP LEV,
    SEXP ordered,
    SEXP teststat,
    SEXP minbucket,
    SEXP lower,
    SEXP give_log
);


/* R_Utils.c */
extern SEXP R_kronecker
(
    SEXP A,
    SEXP B
);
