
/* sum(weights) */
int C_sum_weights
(
    const int *weights,
    const int N
);

/* sum(weights[subset]) */
int C_sum_weights_subset
(
    const int *weights,
    const int N,
    const int *subset,
    const int Nsubset
);

/* colSums(x) */
extern void C_colSums_
(
    const double *x,
    const int N,
    const int P,
    double *P_ans
);

/* rowSums(x) */
extern void C_rowSums_
(
    const double *x,
    const int N,
    const int P,
    double *N_ans
);

/* colSums(x) integer version */
extern void C_colSums_i
(
    const int *x,
    const int N,
    const int P,
    int *P_ans
);

/* rowSums(x) integer version */
extern void C_rowSums_i
(
    const int *x,
    const int N,
    const int P,
    int *N_ans
);

/* colSums(x * weights) */
extern void C_colSums_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
);

/* colSums(x[subsetx,]) */
extern void C_colSums_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
);

/* colSums(x[subsetx,] * weights[subsetx]) */
extern void C_colSums_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
);

/* colSums(x^2) */
extern void C_colSums2_
(
    const double *x,
    const int N,
    const int P,
    double *P_ans
);

/* colSums(x^2 * weights) */
extern void C_colSums2_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
);

/* colSums(x[subsetx,]^2) */
extern void C_colSums2_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
);

/* colSums(x[subsetx,]^2 * weights[subsetx]) */
extern void C_colSums2_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
);

/* colSums((x-center)^2) */
extern void C_colSums2_center_
(
    const double *x,
    const int N,
    const int P,
    const double *centerx,
    double *P_ans
);

/* colSums((x-center)^2 * weights) */
extern void C_colSums2_center_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const double *centerx,
    double *P_ans
);

/* colSums((x[subsetx,] - center)^2) */
extern void C_colSums2_center_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *P_ans
);

/* colSums((x[subsetx,]-center)^2 * weights[subsetx]) */
extern void C_colSums2_center_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *P_ans
);

/* sum_i (t(x[i,]) %*% y[i,]) */
extern void C_KronSums_
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    double *PQ_ans
);

/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
extern void C_KronSums_weights
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    double *PQ_ans
);

/* sum_i (t(x[subsetx[i],]) %*% y[subsety[i],]) */
extern void C_KronSums_subset
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *subsetx,
    const int *subsety,
    const int Nsubset,
    double *PQ_ans
);

/* sum_i weights[subset[i]] (t(x[subset[i],]) %*% y[subset[i],]) */
extern void C_KronSums_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PQ_ans
);

/* sum_i,j weights2d[i, j] * t(x[i,]) %*% y[j,]) */
extern void C_KronSums_2dweights
(
    const double *x,
    const int Lx,
    const int P,
    const double *y,
    const int Ly,
    const int Q,
    const int *weights2d,
    double *PQ_ans
);

/* sum_i (t(x[i,]) %*% x[i,]) */
extern void C_KronSums_sym_
(
    const double *x,
    const int N,
    const int P,
    double *PP_sym_ans
);

/* sum_i weights[i] * (t(x[i,]) %*% x[i,]) */
extern void C_KronSums_sym_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *PP_sym_ans
);

/* sum_i (t(x[subset[i],]) %*% x[subset[i],]) */
extern void C_KronSums_sym_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *PP_sym_ans
);

/* sum_i weights[subset[i]] (t(x[subset[i],]) %*% y[subset[i],]) */
extern void C_KronSums_sym_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PP_sym_ans
);

/* sum_i (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
extern void C_KronSums_sym_center_
(
    const double *x,
    const int N,
    const int P,
    const double *centerx,
    double *PP_sym_ans
);

/* sum_i weights[i] (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
extern void C_KronSums_sym_center_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const double *centerx,
    double *PP_sym_ans
);

/* sum_i (t(x[subset[i],] - centerx) %*% (x[subset[i],] - centerx)) */
extern void C_KronSums_sym_center_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subset,
    const int Nsubset,
    const double *centerx,
    double *PP_sym_ans
);

/* sum_i weights[subsetx[i]] (t(x[subsetx[i],] - centerx) %*%
                            (x[subsetx[i],] - centerx)) */
extern void C_KronSums_sym_center_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *PP_sym_ans
);

/* tapply(1:nrow(y), ix, function(i) colSums(y[i,])) */
extern void C_tapplySum_
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    double *LxQ_ans
);

/* tapply(1:nrow(y), ix, function(i) colSums(weights[i] * y[i,])) */
extern void C_tapplySum_weights
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    const int *weights,
    double *LxQ_ans
);

/* tapply((1:nrow(y))[subsety], ix[subsetx],
          function(i) colSums(y[i,])) */
extern void C_tapplySum_subset
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    const int *subsetx,
    const int *subsety,
    const int Nsubset,
    double *LxQ_ans
);

/* tapply((1:nrow(y))[subset], ix[subset],
          function(i) colSums(weights[i] * y[i,])) */
extern void C_tapplySum_weights_subset
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *LxQ_ans
);

/* sum(weights[i, ] * y[, q]) forall i = 1, ... Lx and q = 0, ..., Q */
extern void C_tapplySum_2d
(
    const double *y,
    const int Ly,
    const int Q,
    const int Lx,
    const int *weights2d,
    double *Lx1Q_ans
);
