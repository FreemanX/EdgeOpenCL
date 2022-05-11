//
// Created by pfxu on 3/15/19.
//

#include "heteroCompLib.h"

void cpuKernel(io_pack *input, io_pack *output, int offset, int length) {
    flushed_printf("cpuKernel executing... ");
    double time = getCurrentTime();
    long a = 1;
    for (long i = 0; i < 10000000; ++i) {
        gettid();
        a += ((1 << 2) + i) / a;
        a += static_cast<long>(exp(i) * pow(exp(i), 12) / log(a));
    }
    flushed_printf("Execution time: %fs \n", getCurrentTime() - time);
}

int main(int argc, char **argv) {
    const int numThread = 8;
    pthread_t thread[numThread];
    cpu_task_pack cpuTaskPack[numThread];
    //size_t dataLength = 10240 * 10240;
    //flushed_printf("Warming up... ");
    //auto *tmp_a = (float *) malloc(dataLength * sizeof(float));
    //auto *tmp_b = (float *) malloc(dataLength * sizeof(float));
    //auto *tmp_c = (float *) malloc(dataLength * sizeof(float));
    //for (int i = 0; i < 100; ++i) {
//  //      vecAddCPU(tmp_a, tmp_b, tmp_c, dataLength);
    //}
    //free(tmp_a);
    //free(tmp_b);
    //free(tmp_c);
    //flushed_printf(" Done!\n");

    for (int i = 0; i < numThread; ++i) {
        task_pack taskPack = {};
        taskPack.taskType = CPU_TASK;
        taskPack.length = 0;
        taskPack.offset = 0;
        cpuTaskPack[i].thread_id = i;
        cpuTaskPack[i].taskPack = &taskPack;
        cpuTaskPack[i].cpuKernel = &cpuKernel;
        cpuTaskPack[i].cpu_id = 7 - i;
        executeOnCPU(&thread[i], &cpuTaskPack[i]);
    }

    for (int i = 0; i < numThread; ++i) {
        pthread_join(thread[i], nullptr);
    }

    //flushed_printf("Second run... \n");

    //for (int i = 0; i < numThread; ++i) {
    //    task_pack taskPack = {};
    //    taskPack.taskType = CPU_TASK;
    //    taskPack.length = 0;
    //    taskPack.offset = 0;
    //    cpuTaskPack[i].thread_id = i;
    //    cpuTaskPack[i].taskPack = &taskPack;
    //    cpuTaskPack[i].cpuKernel = &cpuKernel;
    //    cpuTaskPack[i].cpu_id = 7 - i;
    //    executeOnCPU(&thread[i], &cpuTaskPack[i]);
    //    pthread_join(thread[i], nullptr);
    //}

    return 0;
}

