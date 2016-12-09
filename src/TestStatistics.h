
extern double C_quadform
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *MPinv
);

extern double C_maxtype
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol,
    const int alternative
);

extern void C_standardise
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol
);
