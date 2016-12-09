
extern int NLEVELS
(
    SEXP x
);

extern int NROW
(
    SEXP x
);

extern int NCOL
(
    SEXP x
);

extern void C_kronecker
(
    const double *A,
    const int m,
    const int n,
    const double *B,
    const int r,
    const int s,
    const int overwrite,
    double *ans
);

extern void C_kronecker_sym
(
    const double *A,
    const int m,
    const double *B,
    const int r,
    const int overwrite,
    double *ans
);

extern void C_MPinv_sym
(
    const double *x,
    const int n,
    const double tol,
    double *dMP,
    int *rank
);

extern void rcont2
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
