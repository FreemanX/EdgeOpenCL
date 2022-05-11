//
// Created by pfxu on 4/29/19.
//

//Prototyping hetero-compute concurrent mode

#include "heteroCompLib.h"
#include <algorithm>

using namespace std;

const char *kernelSource = "\n" \
"__kernel void vecAdd(  __global float *a,                       \n" \
"                       __global float *b,                       \n" \
"                       __global float *c,                       \n" \
"                       const int n)                             \n" \
"{                                                               \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    if (id < n)                                                 \n" \
"        c[id] = a[id] + b[id];                                  \n" \
"}                                                               \n";

cl_wrapper *cl;
cl_program program;
cl_kernel kernel;
cl_command_queue exeCmdQueue;
const size_t LOOP_TIME = 4000;


void initCLEnv() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices->deviceId);
    program = cl->createProgram(&kernelSource, 1);
    kernel = cl->createKernel("vecAdd", program);
    exeCmdQueue = cl->createProfilingCmdQueue();
}

void CPUMatAdd(const float *a, const float *b, float *c, int dataLength) {
    for (int i = 0; i < dataLength; ++i) {
//        c[i] = (a[i] + a[i] * b[i]) / b[i];
        c[i] = a[i] + b[i];
    }
}

struct CPU_task {
    int cpu_id;
    float *a;
    float *b;
    float *c;
    int offset;
    int length;
    double timeWasted = 0; // setting affinity may fail, record affinity setting time
};

void *CPUWorker(void *cpu_task) {
    auto *cpuTask = (CPU_task *) cpu_task;
    if (cpuTask->length > 0) {
        int failTime = 0;
        double timer = getCurrentTime();
        cpuTask->timeWasted = 0;
        while (setCurThreadAffinity(cpuTask->cpu_id) != 0) {
            failTime++;
        } // force setting cpu affinity
        if (failTime > 0) {
            cpuTask->timeWasted = getCurrentTime() - timer;
//            flushed_printf("!!!Setting CPU %d affinity failed %d times, costs %f sec\n",
//                           cpuTask->cpu_id, failTime, cpuTask->timeWasted = getCurrentTime() - timer);
        }

        //checkCurThreadAffinity(cpuTask->cpu_id);
        float *wa = cpuTask->a + cpuTask->offset;
        float *wb = cpuTask->b + cpuTask->offset;
        float *wc = cpuTask->c + cpuTask->offset;
        CPUMatAdd(wa, wb, wc, cpuTask->length);
    }
    pthread_exit(nullptr);
}

struct GPU_task {
    cl_mem *a;
    cl_mem *b;
    cl_mem *c;
    size_t offset;
    size_t length;
    size_t dataLength;
};


void *GPUWorker(void *gpu_task) {
    auto *gpuTask = (GPU_task *) gpu_task;

    if (gpuTask->length > 0) {
        cl_wrapper::checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), gpuTask->a));
        cl_wrapper::checkError(clSetKernelArg(kernel, 1, sizeof(cl_mem), gpuTask->b));
        cl_wrapper::checkError(clSetKernelArg(kernel, 2, sizeof(cl_mem), gpuTask->c));
        cl_wrapper::checkError(clSetKernelArg(kernel, 3, sizeof(int), (int *) &gpuTask->dataLength));

        size_t global_size[] = {gpuTask->length};
        size_t wg_size[] = {1024};
        size_t offset[] = {gpuTask->offset};
        cl_event e;
        cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueue,
                                                      kernel,
                                                      1, offset,
                                                      global_size,
                                                      wg_size, 0, nullptr, &e));
        clWaitForEvents(1, &e);
    }
    pthread_exit(nullptr);
}

double executionLayer(size_t *offsets, size_t *workLength, size_t dataLength,
                      ZeroCopyMem<float> *vecA,
                      ZeroCopyMem<float> *vecB,
                      ZeroCopyMem<float> *vecC) {

    pthread_t cpu_worker[8];
    CPU_task cpuTasks[8];
    //launch CPU workers
    for (int i = 0; i < 8; ++i) {
        cpuTasks[i].cpu_id = i;
        cpuTasks[i].a = vecA->hostPtr;
        cpuTasks[i].b = vecB->hostPtr;
        cpuTasks[i].c = vecC->hostPtr;
        cpuTasks[i].offset = offsets[i];
        cpuTasks[i].length = workLength[i];
        if (pthread_create(&cpu_worker[i], nullptr, &CPUWorker, &cpuTasks[i])) {
            flushed_printf("Pthread %d create failed...\n", cpuTasks[i].cpu_id);
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        }
    }

    //launch GPU worker
    pthread_t gpu_worker;
    GPU_task gpuTask{};
    gpuTask.a = &vecA->deviceBuffer;
    gpuTask.b = &vecB->deviceBuffer;
    gpuTask.c = &vecC->deviceBuffer;
    gpuTask.length = workLength[8];
    gpuTask.offset = offsets[8];
    gpuTask.dataLength = dataLength;

    if (pthread_create(&gpu_worker, nullptr, &GPUWorker, &gpuTask)) {
        flushed_printf("GPU thread %d create failed...\n", 8);
        flushed_printf("Error %d: %s\n", errno, strerror(errno));
    }

    for (pthread_t worker : cpu_worker) {
        pthread_join(worker, nullptr);
    }
    pthread_join(gpu_worker, nullptr);
    double affinityCost = 0;
    for (auto &cpuTask : cpuTasks) {
        affinityCost += cpuTask.timeWasted;
    }
    return affinityCost;
}

void schedulingLayer_experimental(ZeroCopyMem<float> *a,
                                  ZeroCopyMem<float> *b,
                                  ZeroCopyMem<float> *c,
                                  int dataLength) {
    float initSplit[3];
    pthread_t testThread;
    int loop_time = std::max(100, (int) LOOP_TIME / (int) log2(dataLength));
    flushed_printf("loop time: %d\n", loop_time);

    auto *tmp_c = (float *) malloc(dataLength * sizeof(float));
    //ZeroCopyMem<float> tmp_cz = init_zero_copy_region<float>(cl, dataLength);
    //int testLength = 10 * 1024;
    CPU_task cpuTask_test{};
    cpuTask_test.a = a->hostPtr;
    cpuTask_test.b = b->hostPtr;
    cpuTask_test.c = tmp_c;
    //cpuTask_test.length = dataLength > testLength ? testLength : dataLength;
    cpuTask_test.length = dataLength;
    cpuTask_test.offset = 0;
    cpuTask_test.cpu_id = 0;

    initSplit[0] = getCurrentTime();
    for (int i = 0; i < loop_time; ++i) {
        if (pthread_create(&testThread, nullptr, &CPUWorker, &cpuTask_test))
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        pthread_join(testThread, nullptr);
    }
    initSplit[0] = (getCurrentTime() - initSplit[0]) / loop_time;
    flushed_printf("Time %d: %f\n", cpuTask_test.cpu_id, initSplit[0]);

    cpuTask_test.cpu_id = 7;
    initSplit[1] = getCurrentTime();
    for (int i = 0; i < loop_time; ++i) {
        if (pthread_create(&testThread, nullptr, &CPUWorker, &cpuTask_test))
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        pthread_join(testThread, nullptr);
    }
    initSplit[1] = (getCurrentTime() - initSplit[1]) / loop_time;
    flushed_printf("Time %d: %f\n", cpuTask_test.cpu_id, initSplit[1]);


    GPU_task gpuTask_test{};
    gpuTask_test.a = &a->deviceBuffer;
    gpuTask_test.b = &b->deviceBuffer;
    //gpuTask_test.c = &tmp_cz.deviceBuffer;
    gpuTask_test.c = &c->deviceBuffer;
    gpuTask_test.offset = 0;
    gpuTask_test.length = cpuTask_test.length;
    gpuTask_test.dataLength = dataLength;
    initSplit[2] = getCurrentTime();
    for (int i = 0; i < loop_time; ++i) {
        if (pthread_create(&testThread, nullptr, &GPUWorker, &gpuTask_test))
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        pthread_join(testThread, nullptr);
    }
    initSplit[2] = (getCurrentTime() - initSplit[2]) / loop_time;
    flushed_printf("Time GPU: %f, diff: %f\n", initSplit[2],
            //               calVecDiff(tmp_cz.hostPtr, tmp_c, cpuTask_test.length));
                   calVecDiff(c->hostPtr, tmp_c, cpuTask_test.length));

    float fraction = -1.0;
    for (float i : initSplit) fraction += 1 / i;
    fraction += 1.0;
    //flushed_printf("fraction = %f\n", fraction);
    float splitSum = 0.0;
    for (int i = 0; i < 3; ++i) {
        initSplit[i] = (1 / initSplit[i]) / fraction;
        flushed_printf("Split %d = %f \n", i, initSplit[i]);
        splitSum += initSplit[i];
    }
    flushed_printf("Split Sum: %f \n\n", splitSum);

    size_t offsets[9];
    size_t workLength[9];
    size_t headIdx = 0;
    for (int i = 0; i < 2; ++i) {
        int length = dataLength * initSplit[i] / 4;
        for (int j = 0; j < 4; ++j) {
            offsets[j + 4 * i] = headIdx;
            workLength[j + 4 * i] = length;
            headIdx += length;
        }
    }
    offsets[8] = headIdx;
    workLength[8] = dataLength - offsets[8];
    flushed_printf("data length: %d, tail: %d\n", dataLength, offsets[8] + workLength[8]);
    for (int i = 0; i < 9; ++i) {
        flushed_printf("head[%d] = %d, workload = %d\n", i, offsets[i], workLength[i]);
    }

    double affCost = 0;
    double timer = getCurrentTime();
    affCost = executionLayer(offsets, workLength, dataLength, a, b, c);
    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer - affCost,
                   calVecDiff(tmp_c, c->hostPtr, dataLength));
    flushed_printf("CPU affinity time: %f sec \n", affCost);

    initSplit[0] = 1;
    initSplit[1] = 0;
    initSplit[2] = 0;
    headIdx = 0;
    for (int i = 0; i < 2; ++i) {
        int length = dataLength * initSplit[i] / 4;
        for (int j = 0; j < 4; ++j) {
            offsets[j + 4 * i] = headIdx;
            workLength[j + 4 * i] = length;
            headIdx += length;
        }
    }
    offsets[8] = headIdx;
    workLength[8] = dataLength - offsets[8];
    flushed_printf("data length: %d, tail: %d\n", dataLength, offsets[8] + workLength[8]);
    for (int i = 0; i < 9; ++i) {
        flushed_printf("head[%d] = %d, workload = %d\n", i, offsets[i], workLength[i]);
    }
    timer = getCurrentTime();
    affCost = executionLayer(offsets, workLength, dataLength, a, b, c);
    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer - affCost,
                   calVecDiff(tmp_c, c->hostPtr, dataLength));
    flushed_printf("CPU affinity time: %f sec \n", affCost);

    initSplit[0] = 0;
    initSplit[1] = 1;
    initSplit[2] = 0;
    headIdx = 0;
    for (int i = 0; i < 2; ++i) {
        int length = dataLength * initSplit[i] / 4;
        for (int j = 0; j < 4; ++j) {
            offsets[j + 4 * i] = headIdx;
            workLength[j + 4 * i] = length;
            headIdx += length;
        }
    }
    offsets[8] = headIdx;
    workLength[8] = dataLength - offsets[8];
    flushed_printf("data length: %d, tail: %d\n", dataLength, offsets[8] + workLength[8]);
    for (int i = 0; i < 9; ++i) {
        flushed_printf("head[%d] = %d, workload = %d\n", i, offsets[i], workLength[i]);
    }

    timer = getCurrentTime();
    affCost = executionLayer(offsets, workLength, dataLength, a, b, c);
    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer - affCost,
                   calVecDiff(tmp_c, c->hostPtr, dataLength));
    flushed_printf("CPU affinity time: %f sec \n", affCost);
//    timer = getCurrentTime();
//    executionLayer(offsets, workLength, dataLength, a, b, c);
//    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer,
//                   calVecDiff(tmp_c, c->hostPtr, dataLength));

    initSplit[0] = 0;
    initSplit[1] = 0;
    initSplit[2] = 1;
    headIdx = 0;
    for (int i = 0; i < 2; ++i) {
        int length = dataLength * initSplit[i] / 4;
        for (int j = 0; j < 4; ++j) {
            offsets[j + 4 * i] = headIdx;
            workLength[j + 4 * i] = length;
            headIdx += length;
        }
    }
    offsets[8] = headIdx;
    workLength[8] = dataLength - offsets[8];
    flushed_printf("data length: %d, tail: %d\n", dataLength, offsets[8] + workLength[8]);
    for (int i = 0; i < 9; ++i) {
        flushed_printf("head[%d] = %d, workload = %d\n", i, offsets[i], workLength[i]);
    }

    timer = getCurrentTime();
    affCost = executionLayer(offsets, workLength, dataLength, a, b, c);
    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer - affCost,
                   calVecDiff(tmp_c, c->hostPtr, dataLength));
    flushed_printf("CPU affinity time: %f sec \n", affCost);
//    timer = getCurrentTime();
//    executionLayer(offsets, workLength, dataLength, a, b, c);
//    flushed_printf("CPU_GPU time = %f, diff: %f \n", getCurrentTime() - timer,
//                   calVecDiff(tmp_c, c->hostPtr, dataLength));
}

int main(int argc, char **argv) {
    const size_t w = 1024;
    const size_t h = 1024;
    const size_t dataLength = w * h * 512 / sizeof(float);
    initCLEnv();
    ZeroCopyMem<float> a = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> b = init_zero_copy_region<float>(cl, dataLength);
    randArrayGenerator<float>(0.0f, 1.0, a.hostPtr, dataLength);
    randArrayGenerator<float>(0.0f, 1.0, b.hostPtr, dataLength);
    ZeroCopyMem<float> c = init_zero_copy_region<float>(cl, dataLength);

    // warm-up GPU
    //cl_wrapper::checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), &a.deviceBuffer));
    //cl_wrapper::checkError(clSetKernelArg(kernel, 1, sizeof(cl_mem), &b.deviceBuffer));
    //cl_wrapper::checkError(clSetKernelArg(kernel, 2, sizeof(cl_mem), &c.deviceBuffer));
    //cl_wrapper::checkError(clSetKernelArg(kernel, 3, sizeof(int), &dataLength));
    //size_t global_size[] = {4096};
    //size_t wg_size[] = {1024};
    //size_t offset[] = {0};
    //cl_event e;
    //cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueue,
    //                                              kernel,
    //                                              1, offset,
    //                                              global_size,
    //                                              wg_size, 0, nullptr, &e));
    //clWaitForEvents(1, &e);
    //--- warm-up GPU


    flushed_printf("schedulingLayer_experimental start... \n");
//    schedulingLayer_experimental(&a, &b, &c, 1);
//    flushed_printf("\n");
//    schedulingLayer_experimental(&a, &b, &c, dataLength / (1024 * 1024));
//    flushed_printf("=====\n");
//    schedulingLayer_experimental(&a, &b, &c, dataLength / 1024);
//    flushed_printf("=====\n");
    schedulingLayer_experimental(&a, &b, &c, dataLength);
    flushed_printf("=====\n");

    return 0;
}
