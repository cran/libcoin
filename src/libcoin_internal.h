
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */
#include <R_ext/stats_package.h> /* for S_rcont2 */

#define ALTERNATIVE_twosided            1
#define ALTERNATIVE_less                2
#define ALTERNATIVE_greater             3

#define TESTSTAT_maximum                1
#define TESTSTAT_quadratic              2

/* S[i, j] for n x n symmetric matrix in lower packed storage allowing for i < j */
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
#define LE(x, y, tol)  ((x) < (y)) || (fabs((x) - (y)) < (tol))
#define GE(x, y, tol)  ((x) > (y)) || (fabs((x) - (y)) < (tol))

#define LinearStatistic_SLOT            0
#define Expectation_SLOT                1
#define Covariance_SLOT                 2
#define Variance_SLOT                   3
#define MPinv_SLOT                      4
#define ExpectationX_SLOT               5
#define varonly_SLOT                    6
#define dim_SLOT                        7
#define ExpectationInfluence_SLOT       8
#define CovarianceInfluence_SLOT        9
#define VarianceInfluence_SLOT          10
#define Xfactor_SLOT                    11
#define Work_SLOT                       12
#define tol_SLOT                        13
#define PermutedLinearStatistic_SLOT    14
#define TableBlock_SLOT                 15
#define Sumweights_SLOT                 16
#define Table_SLOT                      17
