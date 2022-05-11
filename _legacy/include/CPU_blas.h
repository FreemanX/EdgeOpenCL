//
// Created by pfxu on 2/26/19.
//

#ifndef CPU_BLAS_H
#define CPU_BLAS_H

#include "heteroCompLib.h"

double vecAddCPU(const float *a, const float *b, float *c, size_t n);

void gemm_bin(int M, int N, int K, float ALPHA,
              char *A, int lda,
              float *B, int ldb,
              float *C, int ldc);

double gemm_cpu(int TA, int TB, int M, int N, int K, float ALPHA,
                float *A, int lda,
                float *B, int ldb,
                float BETA,
                float *C, int ldc);

#endif
