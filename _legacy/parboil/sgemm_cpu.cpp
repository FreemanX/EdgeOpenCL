//
// Created by pfxu on 7/20/19.
//

#include "heteroCompLib.h"
#include<fstream>
#include<iostream>
#include<vector>

bool readColMajorMatrixFile(const char *fn, int &nr_row, int &nr_col, std::vector<float> &v) {
    std::cerr << "Opening file:" << fn << std::endl;
    std::fstream f(fn, std::fstream::in);
    if (!f.good()) {
        return false;
    }

    // Read # of rows and cols
    f >> nr_row;
    f >> nr_col;

    float data;
    std::cerr << "Matrix dimension: " << nr_row << "x" << nr_col << std::endl;
    while (f.good()) {
        f >> data;
        v.push_back(data);
    }
    v.pop_back(); // remove the duplicated last element
    return true;
}

void
basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb,
           float beta, float *C, int ldc) {
    if ((transa != 'N') && (transa != 'n')) {
        std::cerr << "unsupported value of 'transa' in regtileSgemm()" << std::endl;
        return;
    }

    if ((transb != 'T') && (transb != 't')) {
        std::cerr << "unsupported value of 'transb' in regtileSgemm()" << std::endl;
        return;
    }
#pragma omp parallel for collapse (2)
    for (int mm = 0; mm < m; ++mm) {
        for (int nn = 0; nn < n; ++nn) {
            float c = 0.0f;
            for (int i = 0; i < k; ++i) {
                float a = A[mm + i * lda];
                float b = B[nn + i * ldb];
                c += a * b;
            }
            C[mm + nn * ldc] = C[mm + nn * ldc] * beta + alpha * c;
        }
    }
}

//void runCPU(char *input) {
int runCPU(const char *input, double *exeTime) {
    char matAPath[128];
    char matBPath[128];
    strcpy(matAPath, input);
    strcpy(matBPath, input);
    strcat(matAPath, "matrix1.txt");
    strcat(matBPath, "matrix2.txt");
    int matArow, matAcol;
    int matBrow, matBcol;
    std::vector<float> matA, matBT;
    // load A
    readColMajorMatrixFile(matAPath, matArow, matAcol, matA);
    // load B^T
    readColMajorMatrixFile(matBPath, matBcol, matBrow, matBT);
    // allocate space for C
    std::vector<float> matC(matArow * matBcol);

    flushed_printf("\n\tSgemm Running...\n");
    double timer = getCurrentTime();
    basicSgemm('N', 'T', matArow, matBcol, matAcol, 1.0f, &matA.front(), matArow, &matBT.front(), matBcol, 0.0f,
               &matC.front(), matArow);
    timer = getCurrentTime() - timer;
    flushed_printf("\tSgemm CPU done, time: %f sec\n", timer);
    *exeTime = timer;
    std::cout << "GFLOPs = " << 2. * matArow * matBcol * matAcol / timer / 1e9 << std::endl;
    matA.clear();
    matBT.clear();
    matC.clear();
    return 0;
}

int main(int argc, char **argv) {
    char *input = (char*) "datasets/sgemm/medium/input/";

    std::vector<int> littleCPUList;
    std::vector<int> bigCPUList;
    for (int i = 0; i < 4; ++i) {
        littleCPUList.push_back(i);
        bigCPUList.push_back(i + 4);
    }
    char *input_small = (char*) "datasets/sgemm/small/input/";
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
