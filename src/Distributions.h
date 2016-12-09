
extern double C_chisq_pvalue
(
    const double stat,
    const int df,
    const int lower,
    const int give_log
);

extern double C_perm_pvalue
(
    const int greater,
    const int B,
    const int lower,
    const int give_log
);

extern double C_maxtype_pvalue
(
    const double stat,
    const double *Covariance,
    const int n,
    const int alternative,
    const int lower,
    const int give_log,
    int maxpts,
    double releps,
    double abseps,
    const double tol
);

extern void C_Permute
(
    const int *x,
    const int n,
    int *ans
);

extern void C_PermuteBlock
(
    const int *x,
    const int *table,
    const int Ntable,
    int *ans
);

extern void C_doPermuteBlock
(
    const int *subset,
    const int Nsubset,
    const int *table,
    const int Nlevels,
    int *Nsubset_tmp,
    int *perm
);

extern void C_doPermute
(
    const int *subset,
    const int Nsubset,
    int *Nsubset_tmp,
    int *perm
);

extern void C_setup_subset
(
    const int N, int *N_ans
);

extern void C_setup_subset_weights
(
    const int N,
    const int *weights,
    int *sw_ans
);

extern void C_setup_subset_weights_subset
(
    const int Nsubset,
    const int *weights,
    const int *subset,
    const int *sw_ans
);

extern void C_order_wrt_block
(
    int *subset,
    const int Nsubset,
    const int *block,
    const int *table,
    const int Nlevels
);
