
#ifdef USE_RCONT2_FROM_R
/* This R _header_ file containts the function defintion for S_rcont2 */
#include <R_ext/stats_package.h>
#else
/* use the copy in rcont2.c */
extern void S_rcont2
(
    int *nrow,
    int *ncol,
    int *nrowt,
    int *ncolt,
    int *ntotal,
    double *fact,
    int *jwork,
    int *matrix
);
#endif
