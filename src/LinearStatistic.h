
extern void RC_LinearStatistic_2d
(
    const SEXP x,
    const int N,
    const int P,
    const double *y,
    const int M,
    const int Q,
    const int *weights2d,
    double *PQ_ans
);

extern void RC_LinearStatistic
(
    const SEXP x,
    const int N,
    const int P,
    const double* y,
    const int Q,
    const int *weights,
    const int *sumweights,
    const int *subset,
    const int *Nsubset,
    const int Lb,
    double *PQ_ans
);

extern void RC_PermutedLinearStatistic
(
    const SEXP x,
    const int N,
    const int P,
    const double* y,
    const int Q,
    const int *perm,
    const int *original,
    const int Nperm,
    double *PQ_ans
);

extern void C_ExpectationInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    double *Q_ans
);

extern void C_ExpectationInfluence_weights_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const int *subset,
    const int Nsubset,
    double *Q_ans
);

extern void C_CovarianceInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const double *ExpInf,
    double *QQ_sym_ans
);

extern void C_CovarianceInfluence_weights_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const int *subset,
    const int Nsubset,
    const double *ExpInf,
    double *QQ_sym_ans
);

extern void C_VarianceInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const double *ExpInf,
    double *Q_ans
);

extern void C_VarianceInfluence_weights_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const int *subset,
    const int Nsubset,
    const double *ExpInf,
    double *Q_ans
);

extern void C_ExpectationCoVarianceInfluence
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int *sumweights,
    const int *subset,
    const int *Nsubset,
    const int Lb,
    const int varonly,
    double *LbQ_ans,
    double *LbQ_var_ans,
    double *LbQQ_sym_ans
);

extern void C_ExpectationX_weights
(
    const double* x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
);

extern void C_CovarianceX_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *PP_sym_ans
);

extern void C_CovarianceX_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PP_sym_ans
);

extern void C_VarianceX_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
);

extern void C_ExpectationLinearStatistic
(
    const int P,
    const int Q,
    const double *ExpInf,
    const double *ExpX,
    const int add,
    double *PQ_ans
);

extern void C_CovarianceLinearStatistic
(
    const int P,
    const int Q,
    const double *CovInf,
    const double *ExpX,
    const double *CovX,
    const int sumweights,
    double *PP_sym_tmp,         /* work vector */
    const int add,
    double *PQPQ_sym_ans
);

extern void C_VarianceLinearStatistic
(
    const int P,
    const int Q,
    const double *VarInf,
    const double *ExpX,
    const double *VarX,
    const int sumweights,
    double *P_tmp,              /* work array */
    const int add,
    double *PQ_ans
);

extern void RC_ExpectationCovarianceLinearStatistic
(
    const SEXP x,
    const int N,
    const int P,
    const int Q,
    const int *weights,
    const int *sumweights,
    const int *subset,
    const int *Nsubset,
    const int Lb,
    const double *ExpXtotal,
    const double *ExpInf,
    const double *CovInf,
    double *work,               /* work vector */
    double *PQ_ans,             /* expectation */
    double *PQPQ_sym_ans        /* covariance */
);

extern void RC_ExpectationVarianceLinearStatistic
(
    const SEXP x,
    const int N,
    const int P,
    const int Q,
    const int *weights,
    const int *sumweights,
    const int *subset,
    const int *Nsubset,
    const int Lb,
    const double *ExpXtotal,
    const double *ExpInf,
    const double *VarInf,
    double *work,               /* work vector */
    double *PQ_ans_Exp,         /* expectation */
    double *PQ_ans_Var          /* variance */
);
