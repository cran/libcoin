
#include "libcoin_internal.h"

/* Variables
   x:           a double or integer N x P or Lx x P matrix
                the first row is == 0 in the latter case
   y:           a double N x Q or Ly x Q matrix
                the first two is == 0 in the latter case
   ix:          an integer vector of length N with elements 0...(Lx - 1)
                0 means NA
   weights:     an integer N vector
   weights2d:   an integer Lx x Ly vector
   subset:      an integer Nsubset vector with elements 0...(N - 1)
   subsetx:     an integer Nsubset vector with elements 0...(N - 1)
   subsety:     an integer Nsubset vector with elements 0...(N - 1)
   centerx:     a double P vector centering the columns of x
   centery:     a double Q vector centering the columns of y
   PQ_ans:      return value, a double P x Q matrix
   P_ans:       return value, a double or integer P vector
   N_ans:       return value, a double or integer N vector
   Q_ans:       return value, a double Q vector
   PP_sym_ans:  return value, a symmetric double P x P matrix in lower packed format
   PQ_ans:      return value, a double P x Q matrix
   LxQ_ans:     return value, a double Lx x Q matrix
   Lx1Q_ans:    return value, a double (Lx - 1) x Q matrix
*/


/* sum(weights) */
int C_sum_weights
(
    const int *weights,
    const int N
) {

    long int ans = 0;
    int ret;

    for (int i = 0; i < N; i++) ans += weights[i];

    /* integer overflow can only happen here */
    if (ans > INT_MAX)
        error("sum of weights is larger than INT_MAX");
    ret = (int) ans;

    return(ret);
}

/* sum(weights[subset]) */
int C_sum_weights_subset
(
    const int *weights,
    const int N,
    const int *subset,
    const int Nsubset
) {

    long int ans = 0;
    int ret;

    for (int i = 0; i < Nsubset; i++) ans += weights[subset[i]];

    /* integer overflow can only happen here */
    if (ans > INT_MAX)
        error("sum of weights is larger than INT_MAX");
    ret = (int) ans;

    return(ret);
}

/* colSums(x) */
void C_colSums_
(
    const double *x,
    const int N,
    const int P,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += x[pN + i];
    }
}

/* rowSums(x) */
void C_rowSums_
(
    const double *x,
    const int N,
    const int P,
    double *N_ans
) {

    for (int i = 0; i < N; i++) {
        N_ans[i] = 0.0;
        for (int p = 0; p < P; p++)
            N_ans[i] += x[p * N + i];
    }
}

/* colSums(x) integer version */
void C_colSums_i
(
    const int *x,
    const int N,
    const int P,
    int *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += x[pN + i];
    }
}

/* rowSums(x) integer version */
void C_rowSums_i
(
    const int *x,
    const int N,
    const int P,
    int *N_ans
) {

    for (int i = 0; i < N; i++) {
        N_ans[i] = 0;
        for (int p = 0; p < P; p++)
            N_ans[i] += x[p * N + i];
    }
}


/* colSums(x * weights) */
void C_colSums_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * x[pN + i];
    }
}

/* colSums(x[subsetx,]) */
void C_colSums_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += x[pN + subsetx[i]];
    }
}

/* colSums(x[subsetx,] * weights[subsetx]) */
void C_colSums_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] * x[pN + subsetx[i]];
    }
}

/* colSums(x^2) */
void C_colSums2_
(
    const double *x,
    const int N,
    const int P,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += pow(x[pN + i], 2);
    }
}

/* colSums(x^2 * weights) */
void C_colSums2_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * pow(x[pN + i], 2);
    }
}

/* colSums(x[subsetx,]^2) */
void C_colSums2_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += pow(x[pN + subsetx[i]], 2);
    }
}

/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] * pow(x[pN + subsetx[i]], 2);
    }
}


/* colSums((x-center)^2) */
void C_colSums2_center_
(
    const double *x,
    const int N,
    const int P,
    const double *centerx,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += pow(x[pN + i] - centerx[p], 2);
    }
}

/* colSums((x-center)^2 * weights) */
void C_colSums2_center_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const double *centerx,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * pow(x[pN + i] - centerx[p], 2);
    }
}

/* colSums((x[subsetx,] - center)^2) */
void C_colSums2_center_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += pow(x[pN + subsetx[i]] - centerx[p], 2);
    }
}

/* colSums((x[subsetx,]-center)^2 * weights[subsetx]) */
void C_colSums2_center_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *P_ans
) {

    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] *
                        pow(x[pN + subsetx[i]] - centerx[p], 2);
    }
}


/* sum_i (t(x[i,]) %*% y[i,]) */
void C_KronSums_
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    double *PQ_ans
) {

    int pN, qP, qN;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            for (int i = 0; i < N; i++)
                 PQ_ans[qP + p] +=  y[qN + i] * x[pN + i];
        }
    }
}

/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_weights
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    double *PQ_ans
) {

    int qP, qN;
    double tmp;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = y[qN + i] * weights[i];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += x[p * N + i] * tmp;
        }
    }
}

/* sum_i (t(x[subsetx[i],]) %*% y[subsety[i],]) */
void C_KronSums_subset
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
) {

    int qP, qN, pN, qPp;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            qPp = qP + p;
            for (int i = 0; i < Nsubset; i++)
                PQ_ans[qPp] += y[qN + subsety[i]] * x[pN + subsetx[i]];
        }
    }
}

/* sum_i weights[subset[i]] (t(x[subset[i],]) %*% y[subset[i],]) */
void C_KronSums_weights_subset
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
) {

    int qP, qN;
    double tmp;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = y[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += x[p * N + subset[i]] * tmp;
        }
    }
}

/* sum_i,j weights2d[i, j] * t(x[i,]) %*% y[j,]) */
void C_KronSums_2dweights
(
    const double *x,
    const int Lx,
    const int P,
    const double *y,
    const int Ly,
    const int Q,
    const int *weights2d,
    double *PQ_ans
) {

    int qPp, qLy, pLxi;

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++) {
            PQ_ans[q * P + p] = 0.0;
            qPp = q * P + p;
            qLy = q * Ly;
            for (int i = 0; i < Lx; i++) {
                pLxi = p * Lx + i;
                for (int j = 0; j < Ly; j++)
                      PQ_ans[qPp] += y[qLy + j] * x[pLxi] *
                                     weights2d[j * Lx + i];
            }
        }
    }
}

/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_
(
    const double *x,
    const int N,
    const int P,
    double *PP_sym_ans
) {

    int pN, qN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] +=  x[qN + i] * x[pN + i];
        }
    }
}

/* sum_i weights[i] * (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *PP_sym_ans
) {

    int qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = x[qN + i] * weights[i];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += x[p * N + i] * tmp;
        }
    }
}

/* sum_i (t(x[subset[i],]) %*% x[subset[i],]) */
void C_KronSums_sym_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subsetx,
    const int Nsubset,
    double *PP_sym_ans
) {

    int qN, pN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < Nsubset; i++)
                PP_sym_ans[SpqP] += x[qN + subsetx[i]] * x[pN + subsetx[i]];
        }
    }
}

/* sum_i weights[subset[i]] (t(x[subset[i],]) %*% y[subset[i],]) */
void C_KronSums_sym_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PP_sym_ans
){

    int qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = x[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += x[p * N + subset[i]] * tmp;
        }
    }
}

/* sum_i (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
void C_KronSums_sym_center_
(
    const double *x,
    const int N,
    const int P,
    const double *centerx,
    double *PP_sym_ans
) {

    int qN, pN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] += (x[qN + i] - centerx[q]) *
                                     (x[pN + i] - centerx[p]);
        }
    }
}

/* sum_i weights[i] (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
void C_KronSums_sym_center_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const double *centerx,
    double *PP_sym_ans
) {

    int qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = (x[qN + i] - centerx[q]) * weights[i];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += (x[p * N + i] - centerx[p]) * tmp;
        }
    }
}

/* sum_i (t(x[subset[i],] - centerx) %*% (x[subset[i],] - centerx)) */
void C_KronSums_sym_center_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subset,
    const int Nsubset,
    const double *centerx,
    double *PP_sym_ans
) {

    int qN, pN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < Nsubset; i++)
                PP_sym_ans[SpqP] += (x[qN + subset[i]] - centerx[q]) *
                                    (x[pN + subset[i]] - centerx[p]);

        }
    }
}

/* sum_i weights[subsetx[i]] (t(x[subsetx[i],] - centerx) %*%
                            (x[subsetx[i],] - centerx)) */
void C_KronSums_sym_center_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subsetx,
    const int Nsubset,
    const double *centerx,
    double *PP_sym_ans
) {
    int qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = (x[qN + subsetx[i]] - centerx[q]) * weights[subsetx[i]];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += (x[p * N + subsetx[i]] - centerx[p])
                                           * tmp;
        }
    }
}

/* tapply(1:nrow(y), ix, function(i) colSums(y[i,])) */
void C_tapplySum_
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    double *LxQ_ans
) {

    int qN, qLx, ixi;

    for (int q = 0; q < Lx * Q; q++) LxQ_ans[q] = 0.0;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qLx = q * Lx;
        for (int i = 0; i < N; i++) {
            ixi = ix[i] - 1; /* ix[i] == 0 means NA */
            if (ixi >= 0)
                LxQ_ans[qLx + ixi] += y[qN + i];
        }
    }
}

/* tapply(1:nrow(y), ix, function(i) colSums(weights[i] * y[i,])) */
void C_tapplySum_weights
(
    const double *y,
    const int N,
    const int Q,
    const int *ix,
    const int Lx,
    const int *weights,
    double *LxQ_ans
) {

    int qN, qLx, ixi;

    for (int q = 0; q < Lx * Q; q++) LxQ_ans[q] = 0.0;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qLx = q * Lx;
        for (int i = 0; i < N; i++) {
            ixi = ix[i] - 1;
            if (ixi >= 0)
                LxQ_ans[qLx + ixi] += weights[i] * y[qN + i];
        }
    }
}

/* tapply((1:nrow(y))[subsety], ix[subsetx],
          function(i) colSums(y[i,])) */
void C_tapplySum_subset
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
) {

    int qN, qLx, ixi;

    for (int q = 0; q < Lx * Q; q++) LxQ_ans[q] = 0.0;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qLx = q * Lx;
        for (int i = 0; i < Nsubset; i++) {
            ixi = ix[subsetx[i]] - 1;
            if (ixi >= 0)
                LxQ_ans[qLx + ixi] += y[qN + subsety[i]];
        }
    }
}

/* tapply((1:nrow(y))[subset], ix[subset],
          function(i) colSums(weights[i] * y[i,])) */
void C_tapplySum_weights_subset
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
) {

    int qN, qLx, ixi;

    for (int q = 0; q < Lx * Q; q++) LxQ_ans[q] = 0.0;

    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qLx = q * Lx;
        for (int i = 0; i < Nsubset; i++) {
            ixi = ix[subset[i]] - 1;
            if (ixi >= 0)
                LxQ_ans[qLx + ixi] += weights[subset[i]] * y[qN + subset[i]];
        }
    }
}

/* sum(weights[i, ] * y[, q]) forall i = 1, ... Lx and q = 0, ..., Q */
void C_tapplySum_2d
(
    const double *y,
    const int Ly,
    const int Q,
    const int Lx,
    const int *weights2d,
    double *Lx1Q_ans
) {

    for (int q = 0; q < (Lx - 1) * Q; q++) Lx1Q_ans[q] = 0.0;

    for (int j = 1; j < Ly; j++) { /* j = 0 means NA */
        for (int i = 1; i < Lx; i++) { /* i = 0 means NA */
            for (int q = 0; q < Q; q++)
                Lx1Q_ans[q * (Lx - 1) + (i - 1)] +=
                    weights2d[j * Lx + i] * y[q * Ly + j];
        }
    }
}
