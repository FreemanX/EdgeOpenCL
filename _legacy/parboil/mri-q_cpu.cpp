//
// Created by pfxu on 7/20/19.
//

#include "benchmark.h"
#include <iostream>
#include <endian.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

#if __BYTE_ORDER != __LITTLE_ENDIAN
# error "File I/O is not implemented for this system: wrong endianness."
#endif

extern "C"
void inputData(const char *fName, int *_numK, int *_numX,
               float **kx, float **ky, float **kz,
               float **x, float **y, float **z,
               float **phiR, float **phiI) {
    int numK, numX;
    FILE *fid = fopen(fName, "r");

    if (fid == nullptr) {
        fprintf(stderr, "Cannot open input file\n");
        exit(-1);
    }
    fread(&numK, sizeof(int), 1, fid);
    *_numK = numK;
    fread(&numX, sizeof(int), 1, fid);
    *_numX = numX;
    *kx = (float *) memalign(16, numK * sizeof(float));
    fread(*kx, sizeof(float), numK, fid);
    *ky = (float *) memalign(16, numK * sizeof(float));
    fread(*ky, sizeof(float), numK, fid);
    *kz = (float *) memalign(16, numK * sizeof(float));
    fread(*kz, sizeof(float), numK, fid);
    *x = (float *) memalign(16, numX * sizeof(float));
    fread(*x, sizeof(float), numX, fid);
    *y = (float *) memalign(16, numX * sizeof(float));
    fread(*y, sizeof(float), numX, fid);
    *z = (float *) memalign(16, numX * sizeof(float));
    fread(*z, sizeof(float), numX, fid);
    *phiR = (float *) memalign(16, numK * sizeof(float));
    fread(*phiR, sizeof(float), numK, fid);
    *phiI = (float *) memalign(16, numK * sizeof(float));
    fread(*phiI, sizeof(float), numK, fid);
    fclose(fid);
}

#define PI   3.1415926535897932384626433832795029f
#define PIx2 6.2831853071795864769252867665590058f

struct kValues {
    float Kx;
    float Ky;
    float Kz;
    float PhiMag;
};

inline
void
ComputePhiMagCPU(int numK,
                 float *phiR, float *phiI, float *phiMag) {
    int indexK = 0;
#pragma omp parallel for
    for (indexK = 0; indexK < numK; indexK++) {
        float real = phiR[indexK];
        float imag = phiI[indexK];
        phiMag[indexK] = real * real + imag * imag;
    }
}

inline
void
ComputeQCPU(int numK, int numX,
            struct kValues *kVals,
            float *x, float *y, float *z,
            float *Qr, float *Qi) {
    float expArg;
    float cosArg;
    float sinArg;

    int indexK, indexX;
#pragma omp parallel for
    for (indexK = 0; indexK < numK; indexK++) {
        for (indexX = 0; indexX < numX; indexX++) {
            expArg = PIx2 * (kVals[indexK].Kx * x[indexX] +
                             kVals[indexK].Ky * y[indexX] +
                             kVals[indexK].Kz * z[indexX]);

            cosArg = cosf(expArg);
            sinArg = sinf(expArg);

            float phi = kVals[indexK].PhiMag;
            Qr[indexX] += phi * cosArg;
            Qi[indexX] += phi * sinArg;
        }
    }
}

void createDataStructsCPU(int numK, int numX, float **phiMag,
                          float **Qr, float **Qi) {
    *phiMag = (float *) memalign(16, numK * sizeof(float));
    *Qr = (float *) memalign(16, numX * sizeof(float));
    memset((void *) *Qr, 0, numX * sizeof(float));
    *Qi = (float *) memalign(16, numX * sizeof(float));
    memset((void *) *Qi, 0, numX * sizeof(float));
}

int runCPU(const char *input, double *exeTime) {
    int numX, numK;        /* Number of X and K values */
    int original_numK;        /* Number of K values in input file */
    float *kx, *ky, *kz;        /* K trajectory (3D vectors) */
    float *x, *y, *z;        /* X coordinates (3D vectors) */
    float *phiR, *phiI;        /* Phi values (complex) */
    float *phiMag;        /* Magnitude of Phi */
    float *Qr, *Qi;        /* Q signal (complex) */
    struct kValues *kVals;

    inputData(input, &original_numK, &numX, &kx, &ky, &kz, &x, &y, &z, &phiR, &phiI);
    numK = original_numK;
    printf("%d pixels in output; %d samples in trajectory; using %d samples\n", numX, original_numK, numK);


    flushed_printf("\n\tmri-q Running...\n");
    double timer = getCurrentTime();
    /* Create CPU data structures */
    createDataStructsCPU(numK, numX, &phiMag, &Qr, &Qi);

    ComputePhiMagCPU(numK, phiR, phiI, phiMag);

    kVals = (struct kValues *) calloc(numK, sizeof(struct kValues));
    int k;
#pragma omp parallel for
    for (k = 0; k < numK; k++) {
        kVals[k].Kx = kx[k];
        kVals[k].Ky = ky[k];
        kVals[k].Kz = kz[k];
        kVals[k].PhiMag = phiMag[k];
    }
    ComputeQCPU(numK, numX, kVals, x, y, z, Qr, Qi);
    timer = getCurrentTime() - timer;
    flushed_printf("\tmri-q CPU done, time: %f sec\n", timer);
    *exeTime = timer;
    free(kx);
    free(ky);
    free(kz);
    free(x);
    free(y);
    free(z);
    free(phiR);
    free(phiI);
    free(phiMag);
    free(kVals);
    free(Qr);
    free(Qi);

    return 0;
}

int main(int argc, char **argv) {
    char const *input = "datasets/mri-q/large/input/64_64_64_dataset.bin";
    benchmark(input, nullptr, runCPU);
    return 0;
}


















