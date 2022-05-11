//
// Created by pfxu on 20/8/2019.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_BENCHMARK_H
#define MOBILEHETEROGENOUSPROJECT_CLION_BENCHMARK_H

#include "HeteroComputeAgent.h"

#define OCL_ERRCK_RETVAL cl_wrapper::checkError

void benchmark(const char *input, const char *input_small = nullptr,
               int (*runCPU)(const char *input, double *exeTime) = nullptr,
               int (*runGPU)(const char *input, double *exeTime) = nullptr) {
    std::vector<int> littleCPUList;
    std::vector<int> bigCPUList;
    for (int i = 0; i < 4; ++i) {
        littleCPUList.push_back(i);
        bigCPUList.push_back(i + 4);
    }

    double exetimer[10];
    if (runCPU != nullptr) {
        printf("%s", input);
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

        if (input_small != nullptr) {
            flushed_printf("1 Big: \n");
            setCurThreadAffinity(7);
            omp_set_num_threads(1);
            runCPU(input_small, &exetimer[5]);

            flushed_printf("1 Little: \n");
            setCurThreadAffinity(3);
            omp_set_num_threads(1);
            runCPU(input_small, &exetimer[6]);

            flushed_printf("4 Big: \n");
            omp_set_num_threads(4);
            setCurThreadAffinity(bigCPUList);
            runCPU(input_small, &exetimer[7]);

            flushed_printf("4 Little: \n");
            setCurThreadAffinity(littleCPUList);
            runCPU(input_small, &exetimer[8]);
        }
    }

    if (runGPU != nullptr) {
        flushed_printf("GPU: \n");
        setCurThreadAffinity(bigCPUList);
        runGPU(input, &exetimer[4]);

        if (input_small != nullptr) {
            flushed_printf("GPU: \n");
            setCurThreadAffinity(bigCPUList);
            runGPU(input_small, &exetimer[9]);
        }
    }

    printf("\n Raw data: ");
    for (double i : exetimer) {
        printf("%lf, ", i);
    }
    printf("\n");

    printf("GPU Ratio: \n");
    for (int j = 0; j < 10; ++j) {
        if (j < 5) printf("%lf, ", exetimer[0] / exetimer[j]);
        else if (input_small != nullptr) printf("%lf, ", exetimer[5] / exetimer[j]);
        else printf("%lf, ", 0.0);
    }
    printf("\n");
}

#endif //MOBILEHETEROGENOUSPROJECT_CLION_BENCHMARK_H
