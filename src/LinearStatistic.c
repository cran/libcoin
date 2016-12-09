
#include "libcoin_internal.h"
#include "Sums.h"
#include "Utils.h"
#include "Tables.h"

/* Variables
   x:             a double N x P matrix
   y:             a double N x Q matrix / a double M x Q matrix
   ix:            integer vector of length N with elements 0...(Lx - 1)
   iy:            integer vector of length N with elements 0...(Ly - 1)
   weights:       an integer vector of length N
   subset:        an integer Nsubset vector with elements 0...(N - 1)
   Lb:            number of levels of block
   PQ_ans:        return value, a double P x Q matrix
   P_ans:         return value, a double P vector
   PP_sym_ans:    return value, a symmetric double P x P matrix in lower packed format
   Q_ans:         return value, a double Q vector
   QQ_sym_ans:    return value, a symmetric double Q x Q matrix in lower packed format
   LbQQ_sym_ans:  return values, Lb symmetric double Q x Q matrices in lower packed format
   PQPQ_sym_ans:  return value, a symmetric double PQ x PQ matrix in lower packed format
   ExpInf:        expectation of influence function y, a double Q vector
   CovInf:        covariance of influence function y, a double Q x (Q + 1) / 2 matrix
   ExpX:          "expectation" of x, a double P vector
   CovX:          "covariance" of x, a double P * (P + 1) / 2 matrix
   VarX:          "variance" of x, a double P vector
   PP_tmp:        temp variable, a symmetric double P x P matrix in lower packed format
   P_tmp:         temp variable, a double P vector
   add:           integer; 0 means init return value with 0 and 1 means add to existing values
*/


void C_LinearStatistic_
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    double *PQ_ans
) {

     C_KronSums_(x, N, P, y, Q,
                 PQ_ans);
}

void C_LinearStatistic_weights
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    double *PQ_ans
) {

     C_KronSums_weights(x, N, P, y, Q, weights,
                        PQ_ans);
}

void C_LinearStatistic_subset
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *subset,
    const int Nsubset,
    double *PQ_ans
) {

     C_KronSums_subset(x, N, P, y, Q, subset, subset, Nsubset,
                       PQ_ans);
}

void C_LinearStatistic_weights_subset
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

     C_KronSums_weights_subset(x, N, P, y, Q, weights, subset, Nsubset,
                               PQ_ans);
}


void C_PermutedLinearStatistic_
(
    const double *x,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *perm,
    const int *original,
    const int Nperm,
    double *PQ_ans
) {

     C_KronSums_subset(x, N, P, y, Q, perm, original, Nperm,
                       PQ_ans);
}

void RC_LinearStatistic_2d
(
    const SEXP x,
    const int N,
    const int P,
    const double *y,
    const int M,
    const int Q,
    const int *weights2d,
    double *PQ_ans
) {

    if (isInteger(x)) {
        C_tapplySum_2d(y, M, Q, P + 1, weights2d,
                       PQ_ans);
    } else {
        C_KronSums_2dweights(REAL(x), N, P, y, M, Q, weights2d,
                             PQ_ans);
    }
}

void C_LinearStatisticXfactor_
(
    const int *ix,
    const int N,
    const int P,
    const double *y,
    const int Q,
    double *PQ_ans
) {

    C_tapplySum_(y, N, Q, ix, P,
                 PQ_ans);
}

void C_LinearStatisticXfactor_weights
(
    const int *ix,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    double *PQ_ans
) {

    C_tapplySum_weights(y, N, Q, ix, P, weights,
                        PQ_ans);
}

void C_LinearStatisticXfactor_subset
(
    const int *ix,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *subset,
    const int Nsubset,
    double *PQ_ans
) {

    C_tapplySum_subset(y, N, Q, ix, P, subset, subset, Nsubset,
                       PQ_ans);
}

void C_LinearStatisticXfactor_weights_subset
(
    const int *ix,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PQ_ans
) {

    C_tapplySum_weights_subset(y, N, Q, ix, P, weights, subset, Nsubset,
                               PQ_ans);
}

void C_PermutedLinearStatisticXfactor_
(
    const int *ix,
    const int N,
    const int P,
    const double *y,
    const int Q,
    const int *perm,
    const int *original,
    const int Nperm,
    double *PQ_ans
) {

    C_tapplySum_subset(y, N, Q, ix, P, perm, original, Nperm,
                       PQ_ans);
}

void C_LinearStatisticXfactor
(
    const int *x,
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
) {

    int sw = 0, ns = 0;

    for (int b = 0; b < Lb; b++) {
        sw = sw + sumweights[b];
        ns = ns + Nsubset[b];
    }

    if (ns < 0) {  /* means: no subset given */
        if (sw < 0) {  /* means: no weights given */
              C_LinearStatisticXfactor_(x, N, P, y, Q,
                                        PQ_ans);
        } else {
              C_LinearStatisticXfactor_weights(x, N, P, y, Q, weights,
                                               PQ_ans);
        }
    } else {
        if (sw < 0) {
            C_LinearStatisticXfactor_subset(x, N, P, y, Q, subset, ns,
                                            PQ_ans);
        } else {
            C_LinearStatisticXfactor_weights_subset(x, N, P, y, Q, weights,
                                                    subset, ns,
                                                    PQ_ans);
        }
    }
}

void RC_LinearStatistic
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
) {

    int sw = 0, ns = 0;

    if (isInteger(x)) {
        C_LinearStatisticXfactor(INTEGER(x), N, P, y, Q, weights,
                                 sumweights, subset, Nsubset, Lb,
                                 PQ_ans);
    } else {
        for (int b = 0; b < Lb; b++) {
            sw = sw + sumweights[b];
            ns = ns + Nsubset[b];
        }

        if (ns < 0) {  /* means: no subset given */
            if (sw < 0) {  /* means: no weights given */
                  C_LinearStatistic_(REAL(x), N, P, y, Q, PQ_ans);
            } else {
                  C_LinearStatistic_weights(REAL(x), N, P, y, Q, weights,
                                            PQ_ans);
        }
        } else {
            if (sw < 0) {
                C_LinearStatistic_subset(REAL(x), N, P, y, Q, subset, ns,
                                         PQ_ans);
            } else {
                C_LinearStatistic_weights_subset(REAL(x), N, P, y, Q, weights,
                                                 subset, ns,
                                                 PQ_ans);
            }
        }
    }
}

void RC_PermutedLinearStatistic
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
) {

    if (isInteger(x)) {
            C_PermutedLinearStatisticXfactor_(INTEGER(x), N, P, y, Q,
                                              perm, original, Nperm,
                                              PQ_ans);
    } else {
            C_PermutedLinearStatistic_(REAL(x), N, P, y, Q,
                                       perm, original, Nperm,
                                       PQ_ans);
    }
}

void C_ExpectationInfluence_
(
    const double* y,
    const int N,
    const int Q,
    double *Q_ans
) {

     C_colSums_(y, N, Q,
                Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_ExpectationInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    double *Q_ans
) {

     C_colSums_weights(y, N, Q, weights,
                       Q_ans);
     if (sumweights > 0) {
         for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
     }
}

void C_ExpectationInfluence_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *subset,
    const int Nsubset,
    double *Q_ans
) {

     C_colSums_subset(y, N, Q, subset, Nsubset,
                      Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_ExpectationInfluence_weights_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const int *subset,
    const int Nsubset,
    double *Q_ans
) {

     C_colSums_weights_subset(y, N, Q, weights, subset, Nsubset,
                              Q_ans);
     if (sumweights > 0) {
         for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
     }
}


void C_CovarianceInfluence_
(
    const double* y,
    const int N,
    const int Q,
    const double *ExpInf,
    double *QQ_sym_ans
) {

     C_KronSums_sym_center_(y, N, Q, ExpInf,
                            QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++)
         QQ_sym_ans[q] = QQ_sym_ans[q] / N;
}

void C_CovarianceInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const double *ExpInf,
    double *QQ_sym_ans
) {

     C_KronSums_sym_center_weights(y, N, Q, weights, ExpInf,
                                   QQ_sym_ans);
     if (sumweights > 0) {
         for (int q = 0; q < Q * (Q + 1) / 2; q++)
             QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
     }
}

void C_CovarianceInfluence_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *subset,
    const int Nsubset,
    const double *ExpInf,
    double *QQ_sym_ans
) {

     C_KronSums_sym_center_subset(y, N, Q, subset, Nsubset, ExpInf,
                                  QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++)
         QQ_sym_ans[q] = QQ_sym_ans[q] / Nsubset;
}

void C_CovarianceInfluence_weights_subset
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
) {

     C_KronSums_sym_center_weights_subset(y, N, Q, weights, subset, Nsubset,
                                          ExpInf,
                                          QQ_sym_ans);
     if (sumweights > 0) {                                          
         for (int q = 0; q < Q * (Q + 1) / 2; q++)
             QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
     }
}

void C_VarianceInfluence_
(
    const double* y,
    const int N,
    const int Q,
    const double *ExpInf,
    double *Q_ans
) {

     C_colSums2_center_(y, N, Q, ExpInf,
                        Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_VarianceInfluence_weights
(
    const double* y,
    const int N,
    const int Q,
    const int *weights,
    const int sumweights,
    const double *ExpInf,
    double *Q_ans
) {

     C_colSums2_center_weights(y, N, Q, weights, ExpInf,
                               Q_ans);
     if (sumweights > 0) {                               
         for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
     }
}

void C_VarianceInfluence_subset
(
    const double* y,
    const int N,
    const int Q,
    const int *subset,
    const int Nsubset,
    const double *ExpInf,
    double *Q_ans
) {

     C_colSums2_center_subset(y, N, Q, subset, Nsubset, ExpInf,
                              Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_VarianceInfluence_weights_subset
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
) {

     C_colSums2_center_weights_subset(y, N, Q, weights, subset, Nsubset,
                                      ExpInf,
                                      Q_ans);
     if (sumweights > 0) {                                      
         for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
     }
}

void C_ExpectationCoVarianceInfluence
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
) {

     int ns = 0;
     double *ExpInf, *CovInf, *VarInf;

     for (int b = 0; b < Lb; b++) {
         ExpInf = LbQ_ans + b * Q;
         VarInf = LbQ_var_ans + b * Q;
         CovInf = LbQQ_sym_ans + b * Q * (Q + 1) / 2;
         if (Nsubset[b] < 0) {  /* means: no subset given */
             if (sumweights[b] < 0) {  /* means: no weights given */
                 C_ExpectationInfluence_(y, N, Q,
                                         ExpInf);
                 /* compute both for the time being (for later reuse) */
                 C_VarianceInfluence_(y, N, Q, ExpInf,
                                      VarInf);
                 C_CovarianceInfluence_(y, N, Q, ExpInf,
                                        CovInf);
             } else {
                 C_ExpectationInfluence_weights(y, N, Q, weights, sumweights[b],
                                                ExpInf);
                 /* compute both for the time being (for later reuse) */
                 C_VarianceInfluence_weights(y, N, Q, weights,
                                             sumweights[b], ExpInf,
                                             VarInf);
                 C_CovarianceInfluence_weights(y, N, Q, weights,
                                               sumweights[b], ExpInf,
                                               CovInf);
             }
         } else {
             if (sumweights[b] < 0) {
                 C_ExpectationInfluence_subset(y, N, Q, subset + ns, Nsubset[b],
                                               ExpInf);
                 /* compute both for the time being (for later reuse) */
                 C_VarianceInfluence_subset(y, N, Q, subset + ns,
                                            Nsubset[b], ExpInf,
                                            VarInf);
                 C_CovarianceInfluence_subset(y, N, Q, subset + ns,
                                              Nsubset[b], ExpInf,
                                              CovInf);
             } else {
                 C_ExpectationInfluence_weights_subset(y, N, Q, weights,
                                                       sumweights[b], subset + ns,
                                                       Nsubset[b],
                                                       ExpInf);
                 /* compute both for the time being (for later reuse) */
                 C_VarianceInfluence_weights_subset(y, N, Q, weights,
                                                    sumweights[b], subset + ns,
                                                    Nsubset[b], ExpInf,
                                                    VarInf);
                 C_CovarianceInfluence_weights_subset(y, N, Q, weights,
                                                      sumweights[b],
                                                      subset + ns,
                                                      Nsubset[b],
                                                      ExpInf,
                                                      CovInf);
             }
             ns = ns + Nsubset[b];
         }
     }
}


void C_ExpectationX_
(
    const double* x,
    const int N,
    const int P,
    double *P_ans
) {

     C_colSums_(x, N, P,
                P_ans);
}

void C_ExpectationX_weights
(
    const double* x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
) {

     C_colSums_weights(x, N, P, weights,
                       P_ans);
}

void C_ExpectationX_subset
(
    const double* x,
    const int N,
    const int P,
    const int *subset,
    const int Nsubset,
    double *P_ans
) {

     C_colSums_subset(x, N, P, subset, Nsubset,
                      P_ans);
}

void C_ExpectationX_weights_subset
(
    const double* x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *P_ans
) {

     C_colSums_weights_subset(x, N, P, weights, subset, Nsubset,
                              P_ans);
}

void C_CovarianceX_
(
    const double *x,
    const int N,
    const int P,
    double *PP_sym_ans
) {

     C_KronSums_sym_(x, N, P,
                     PP_sym_ans);
}

void C_CovarianceX_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *PP_sym_ans
) {

     C_KronSums_sym_weights(x, N, P, weights,
                            PP_sym_ans);
}

void C_CovarianceX_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subset,
    const int Nsubset,
    double *PP_sym_ans
) {

     C_KronSums_sym_subset(x, N, P, subset, Nsubset,
                           PP_sym_ans);
}

void C_CovarianceX_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *PP_sym_ans
) {

     C_KronSums_sym_weights_subset(x, N, P, weights, subset, Nsubset,
                                   PP_sym_ans);
}

void C_VarianceX_
(
    const double *x,
    const int N,
    const int P,
    double *P_ans
) {

     C_colSums2_(x, N, P,
                 P_ans);
}

void C_VarianceX_weights
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    double *P_ans
) {

     C_colSums2_weights(x, N, P, weights,
                        P_ans);
}

void C_VarianceX_subset
(
    const double *x,
    const int N,
    const int P,
    const int *subset,
    const int Nsubset,
    double *P_ans
) {

     C_colSums2_subset(x, N, P, subset, Nsubset,
                       P_ans);
}

void C_VarianceX_weights_subset
(
    const double *x,
    const int N,
    const int P,
    const int *weights,
    const int *subset,
    const int Nsubset,
    double *P_ans
) {

     C_colSums2_weights_subset(x, N, P, weights, subset, Nsubset,
                               P_ans);
}

void C_ExpectationLinearStatistic
(
    const int P,
    const int Q,
    const double *ExpInf,
    const double *ExpX,
    const int add,
    double *PQ_ans
) {

    if (!add)
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] += ExpX[p] * ExpInf[q];
    }
}

void C_CovarianceLinearStatistic
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
) {

    double f1 = (double) sumweights / (sumweights - 1);
    double f2 = 1.0 / (sumweights - 1);
    double tmp;

    if (P * Q == 1) {
        tmp = f1 * CovInf[0] * CovX[0];
        tmp -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
        if (add) {
            PQPQ_sym_ans[0] += tmp;
        } else {
            PQPQ_sym_ans[0] = tmp;
        }
    } else {
        C_KronSums_sym_(ExpX, 1, P,
                        PP_sym_tmp);
        for (int p = 0; p < P * (P + 1) / 2; p++)
            PP_sym_tmp[p] = f1 * CovX[p] - f2 * PP_sym_tmp[p];
        C_kronecker_sym(CovInf, Q, PP_sym_tmp, P, 1 - (add >= 1),
                        PQPQ_sym_ans);
    }
}

void C_VarianceLinearStatistic
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
) {

    if (P * Q == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX,
                                    sumweights, P_tmp, (add >= 1),
                                    PQ_ans);
    } else {

        double f1 = (double) sumweights / (sumweights - 1);
        double f2 = 1.0 / (sumweights - 1);
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, 1 - (add >= 1),
                    PQ_ans);
    }
}

void RC_ExpectationCovarianceLinearStatistic
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
    double *ExpXtotal,          /* E(g(X)) */
    const double *ExpInf,
    const double *CovInf,
    double *work,               /* work vector */
    double *PQ_ans,             /* expectation */
    double *PQPQ_sym_ans        /* covariance */
) {

     int ns = 0, sw = 0, *ix;
     const int *subtmp;
     double *ExpX, *CovX, *PPtmp;

     /* work[0] counts NAs in ix (ix[i] == 0)  */
     ExpX = work + 1;
     CovX = ExpX + P;
     PPtmp = CovX + P * (P + 1) / 2;

     for (int p = 0; p < P; p++) ExpXtotal[p] = 0.0;

     for (int b = 0; b < Lb; b++) {
         for (int i = 0; i < P + 2 * P * (P + 1) / 2 + 1; i++)
             work[i] = 0.0;
         if (Nsubset[b] < 0) {  /* means: no subset given */
             if (sumweights[b] < 0) {  /* means: no weights given */
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     /* work[0] counts NAs */
                     for (int i = 0; i < N; i++) work[ix[i]]++;
                     /* CovX = diag(ExpX) */
                     for (int p = 0; p < P; p++)
                         CovX[S(p, p, P)] = ExpX[p];
                 } else {
                     C_ExpectationX_(REAL(x), N, P,
                                     ExpX);
                     C_CovarianceX_(REAL(x), N, P,
                                    CovX);
                 }
                 sw = N;
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     for (int i = 0; i < N; i++)
                         work[ix[i]] += (double) weights[i];
                     for (int p = 0; p < P; p++)
                         CovX[S(p, p, P)] = ExpX[p];
                 } else {
                     C_ExpectationX_weights(REAL(x), N, P, weights,
                                            ExpX);
                     C_CovarianceX_weights(REAL(x), N, P, weights,
                                           CovX);
                 }
                 sw = sumweights[b];
             }
         } else {
             if (sumweights[b] < 0) {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++)
                         work[ix[subtmp[i]]]++;
                     for (int p = 0; p < P; p++)
                          CovX[S(p, p, P)] = ExpX[p];
                 } else {
                     C_ExpectationX_subset(REAL(x), N, P, subset + ns,
                                           Nsubset[b],
                                           ExpX);
                     C_CovarianceX_subset(REAL(x), N, P, subset + ns,
                                          Nsubset[b],
                                          CovX);
                 }
                 sw = Nsubset[b];
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++)
                         work[ix[subtmp[i]]] += (double) weights[subtmp[i]];
                     for (int p = 0; p < P; p++)
                          CovX[S(p, p, P)] = ExpX[p];
                 } else {
                     C_ExpectationX_weights_subset(REAL(x), N, P, weights,
                                                   subset + ns, Nsubset[b],
                                                   ExpX);
                     C_CovarianceX_weights_subset(REAL(x), N, P, weights,
                                                  subset + ns, Nsubset[b],
                                                  CovX);
                 }
                 sw = sumweights[b];
             }
         }
         for (int p = 0; p < P; p++) ExpXtotal[p] += ExpX[p];

         C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b,
                                      PQ_ans);
         C_CovarianceLinearStatistic(P, Q, CovInf + b * Q * (Q + 1) / 2,
                                     ExpX, CovX, sw, PPtmp, b,
                                     PQPQ_sym_ans);
         ns = ns + Nsubset[b];
     }
}

void RC_ExpectationVarianceLinearStatistic
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
    double *ExpXtotal,          /* E(g(X)) */
    const double *ExpInf,
    const double *VarInf,
    double *work,               /* work vector */
    double *PQ_ans_Exp,         /* expectation */
    double *PQ_ans_Var          /* variance */
) {

     int ns = 0, sw = 0, *ix;
     const int *subtmp;
     double *ExpX, *VarX, *PPtmp;

     /* work[0] counts NAs in ix (ix[i] = 0)*/
     ExpX = work + 1;
     VarX = ExpX + P;
     PPtmp = VarX + P;

     for (int p = 0; p < P; p++) ExpXtotal[p] = 0.0;

     for (int b = 0; b < Lb; b++) {
         for (int i = 0; i < 3 * P + 1; i++) work[i] = 0.0;
         if (Nsubset[b] < 0) {  /* means: no subset given */
             if (sumweights[b] < 0) {  /* means: no weights given */
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     /* work[0] counts NAs */
                     for (int i = 0; i < N; i++) work[ix[i]]++;
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_(REAL(x), N, P,
                                     ExpX);
                     C_VarianceX_(REAL(x), N, P,
                                  VarX);
                 }
                 sw = N;
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     for (int i = 0; i < N; i++)
                         work[ix[i]] += (double) weights[i];
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_weights(REAL(x), N, P, weights,
                                            ExpX);
                     C_VarianceX_weights(REAL(x), N, P, weights,
                                         VarX);
                 }
                 sw = sumweights[b];
             }
         } else {
             if (sumweights[b] < 0) {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++)
                         work[ix[subtmp[i]]]++;
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_subset(REAL(x), N, P, subset + ns,
                                           Nsubset[b],
                                           ExpX);
                     C_VarianceX_subset(REAL(x), N, P, subset + ns,
                                        Nsubset[b],
                                        VarX);
                 }
                 sw = Nsubset[b];
             } else {
                 if (isInteger(x)) {
                     ix = INTEGER(x);
                     subtmp = subset + ns;
                     for (int i = 0; i < Nsubset[b]; i++)
                         work[ix[subtmp[i]]] += (double) weights[subtmp[i]];
                     VarX = ExpX;
                 } else {
                     C_ExpectationX_weights_subset(REAL(x), N, P, weights,
                                                   subset + ns, Nsubset[b],
                                                   ExpX);
                     C_VarianceX_weights_subset(REAL(x), N, P, weights,
                                                subset + ns, Nsubset[b],
                                                VarX);
                 }
                 sw = sumweights[b];
             }
         }

         for (int p = 0; p < P; p++) ExpXtotal[p] += ExpX[p];

         C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b,
                                      PQ_ans_Exp);
         C_VarianceLinearStatistic(P, Q, VarInf + b * Q, ExpX, VarX,
                                   sw, PPtmp, b,
                                   PQ_ans_Var);
         ns = ns + Nsubset[b];
     }
}
