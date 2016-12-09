
void C_ordered_Xfactor_block
(
    const double *linstat,
    const double *expect,
    const double *covar,
    const int P,
    const int Q,
    const double *ExpX,
    const int B,
    const double* blinstat,
    const int minbucket,
    const double tol,
    const int teststat,
    int *wmax,
    double *maxstat,
    double *pval,
    const int lower,
    const int give_log
);

void C_ordered_Xfactor
(
    const double *linstat,
    const double *expect,
    const double *varinf,
    const double *covinf,
    const int P,
    const int Q,
    const double *ExpX,
    const int B,
    const double *blinstat,
    const int minbucket,
    const double tol,
    const int teststat,
    int *wmax,
    double *maxstat,
    double *pval,
    const int lower,
    const int give_log
);

void C_unordered_Xfactor_block
(
    const double *linstat,
    const double *expect,
    const double *covar,
    const int P,
    const int Q,
    const double *ExpX,
    const int B,
    const double* blinstat,
    const int minbucket,
    const double tol,
    const int teststat,
    int *wmax,
    double *maxstat,
    double *pval,
    const int lower,
    const int give_log
);

void C_unordered_Xfactor
(
    const double *linstat,
    const double *expect,
    const double *varinf,
    const double *covinf,
    const int P,
    const int Q,
    const double *ExpX,
    const int B,
    const double *blinstat,
    const int minbucket,
    const double tol,
    const int teststat,
    int *wmax,
    double *maxstat,
    double *pval,
    const int lower,
    const int give_log
);
