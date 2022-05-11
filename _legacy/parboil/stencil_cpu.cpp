//
// Created by pfxu on 7/20/19.
//
#include "HeteroComputeAgent.h"

#define Index3D(_nx, _ny, _i, _j, _k) ((_i)+_nx*((_j)+_ny*(_k)))

static int read_data(float *A0, int nx, int ny, int nz, FILE *fp) {
    int s = 0;
    int i, j, k;
    for (i = 0; i < nz; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nx; k++) {
                fread(A0 + s, sizeof(float), 1, fp);
                s++;
            }
        }
    }
    return 0;
}

void cpu_stencil(float c0, float c1, float *A0, float *Anext, const int nx, const int ny, const int nz) {
    int i;
#pragma omp parallel for
    for (i = 1; i < nx - 1; i++) {
        int j, k;
        for (j = 1; j < ny - 1; j++) {
            for (k = 1; k < nz - 1; k++) {
                //i      #pragma omp critical
                Anext[Index3D(nx, ny, i, j, k)] =
                        (A0[Index3D(nx, ny, i, j, k + 1)] +
                         A0[Index3D(nx, ny, i, j, k - 1)] +
                         A0[Index3D(nx, ny, i, j + 1, k)] +
                         A0[Index3D(nx, ny, i, j - 1, k)] +
                         A0[Index3D(nx, ny, i + 1, j, k)] +
                         A0[Index3D(nx, ny, i - 1, j, k)]) * c1 - A0[Index3D(nx, ny, i, j, k)] * c0;
            }
        }
    }

}

//void runCPU(char *input) {
int runCPU(const char *input, double *exeTime) {
    int nx, ny, nz;
    int size;
    int iteration;
    float c0 = 1.0f / 6.0f;
    float c1 = 1.0f / 6.0f / 6.0f;

    nx = 128;
    ny = 128;
    nz = 32;
    iteration = 1000;
    //host data
    float *h_A0;
    float *h_Anext;
    size = nx * ny * nz;

    h_A0 = (float *) malloc(sizeof(float) * size);
    h_Anext = (float *) malloc(sizeof(float) * size);
    FILE *fp = fopen(input, "rb");
    read_data(h_A0, nx, ny, nz, fp);
    fclose(fp);
    memcpy(h_Anext, h_A0, sizeof(float) * size);

    flushed_printf("\n\tstencil Running...\n");
    double timer = getCurrentTime();

    int t;
    for (t = 0; t < iteration; t++) {
        cpu_stencil(c0, c1, h_A0, h_Anext, nx, ny, nz);
        float *temp = h_A0;
        h_A0 = h_Anext;
        h_Anext = temp;
    }
    float *temp = h_A0;
    h_A0 = h_Anext;
    h_Anext = temp;

    timer = getCurrentTime() - timer;
    flushed_printf("\tstencil CPU done, time: %f sec\n", timer);
    *exeTime = timer;

    free(h_A0);
    free(h_Anext);
    return 0;
}

int main(int argc, char **argv) {
    char *input = "datasets/stencil/default/input/512x512x64x100.bin";
    std::vector<int> littleCPUList;
    std::vector<int> bigCPUList;
    for (int i = 0; i < 4; ++i) {
        littleCPUList.push_back(i);
        bigCPUList.push_back(i + 4);
    }
    char *input_small = "datasets/stencil/small/input/128x128x32.bin";
    double exetimer[8];

    flushed_printf("1 Big: \n");
    setCurThreadAffinity(7);
    omp_set_num_threads(1);
    runCPU(input, &exetimer[0]);

    flushed_printf("1 Little: \n");
    setCurThreadAffinity(3);
    omp_set_num_threads(1);
    runCPU(input, &exetimer[1]);

    flushed_printf("4 Big: \n");
    omp_set_num_threads(4);
    setCurThreadAffinity(bigCPUList);
    runCPU(input, &exetimer[2]);

    flushed_printf("4 Little: \n");
    setCurThreadAffinity(littleCPUList);
    runCPU(input, &exetimer[3]);


    flushed_printf("1 Big: \n");
    setCurThreadAffinity(7);
    omp_set_num_threads(1);
    runCPU(input_small, &exetimer[4]);

    flushed_printf("1 Little: \n");
    setCurThreadAffinity(3);
    omp_set_num_threads(1);
    runCPU(input_small, &exetimer[5]);

    flushed_printf("4 Big: \n");
    omp_set_num_threads(4);
    setCurThreadAffinity(bigCPUList);
    runCPU(input_small, &exetimer[6]);

    flushed_printf("4 Little: \n");
    setCurThreadAffinity(littleCPUList);
    runCPU(input_small, &exetimer[7]);

    printf("\n Raw data: ");
    for (int i = 0; i < 8; ++i) {
        printf("%lf, ", exetimer[i]) ;
    }
    printf("\n Ratio: ");
    for (int i = 0; i < 8; i+=2) {
        printf("%lf, ", exetimer[i + 1]/exetimer[i]) ;
    }
    printf("\n");
    return 0;
}

















