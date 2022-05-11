//
// Created by pfxu on 2/26/19.
//

#include "CPU_blas.h"

double vecAddCPU(const float *a, const float *b, float *c, const size_t n) {
    int i;
    double timer = getCurrentTime();
#pragma omp parallel for default(none) private(i) shared(a, b, c)
    for (i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
    return getCurrentTime() - timer;
}

void gemm_bin(int M, int N, int K, float ALPHA,
              char *A, int lda,
              float *B, int ldb,
              float *C, int ldc) {
    int i, j, k;
    for (i = 0; i < M; ++i) {
        for (k = 0; k < K; ++k) {
            char A_PART = A[i * lda + k];
            if (A_PART) {
                for (j = 0; j < N; ++j) {
                    C[i * ldc + j] += B[k * ldb + j];
                }
            } else {
                for (j = 0; j < N; ++j) {
                    C[i * ldc + j] -= B[k * ldb + j];
                }
            }
        }
    }
}

void gemm_nn(int M, int N, int K, float ALPHA,
             float *A, int lda,
             float *B, int ldb,
             float *C, int ldc) {
    int i, j, k;
    omp_set_dynamic(0);
#pragma omp parallel for private(k, j)
    for (i = 0; i < M; ++i) {
        for (k = 0; k < K; ++k) {
            register float A_PART = ALPHA * A[i * lda + k];
            for (j = 0; j < N; ++j) {
                C[i * ldc + j] += A_PART * B[k * ldb + j];
            }
        }
    }
}

void gemm_nt(int M, int N, int K, float ALPHA,
             float *A, int lda,
             float *B, int ldb,
             float *C, int ldc) {
    int i, j, k;
#pragma omp parallel for private(k, j)
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            register float sum = 0;
            for (k = 0; k < K; ++k) {
                sum += ALPHA * A[i * lda + k] * B[j * ldb + k];
            }
            C[i * ldc + j] += sum;
        }
    }
}

void gemm_tn(int M, int N, int K, float ALPHA,
             float *A, int lda,
             float *B, int ldb,
             float *C, int ldc) {
    int i, j, k;
//    #pragma omp parallel for
#pragma omp parallel for private(k, j)
    for (i = 0; i < M; ++i) {
        for (k = 0; k < K; ++k) {
            register float A_PART = ALPHA * A[k * lda + i];
            for (j = 0; j < N; ++j) {
                C[i * ldc + j] += A_PART * B[k * ldb + j];
            }
        }
    }
}

void gemm_tt(int M, int N, int K, float ALPHA,
             float *A, int lda,
             float *B, int ldb,
             float *C, int ldc) {
    int i, j, k;
//    #pragma omp parallel for
#pragma omp parallel for private(k, j)
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            register float sum = 0;
            for (k = 0; k < K; ++k) {
                sum += ALPHA * A[i + k * lda] * B[k + j * ldb];
            }
            C[i * ldc + j] += sum;
        }
    }
}

// C = ALPHA*Op(A) * Op(B) + BETA*C, where A(m, k) B(k, n) C(m, n)
double gemm_cpu(int TA, int TB, int M, int N, int K, float ALPHA,
                float *A, int lda,
                float *B, int ldb,
                float BETA,
                float *C, int ldc) {
    double time = getCurrentTime();
    int i, j;
#pragma omp parallel for private(i, j)
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            C[i * ldc + j] *= BETA;
        }
    }
    if (!TA && !TB)
        gemm_nn(M, N, K, ALPHA, A, lda, B, ldb, C, ldc);
    else if (TA && !TB)
        gemm_tn(M, N, K, ALPHA, A, lda, B, ldb, C, ldc);
    else if (!TA && TB)
        gemm_nt(M, N, K, ALPHA, A, lda, B, ldb, C, ldc);
    else
        gemm_tt(M, N, K, ALPHA, A, lda, B, ldb, C, ldc);

    return getCurrentTime() - time;
}

