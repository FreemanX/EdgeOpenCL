//
// Created by pfxu on 3/18/19.
//
#include "heteroCompLib.h"

const char *kernelSource = "\n" \
"__kernel void vecAdd(  __global float *a,                       \n" \
"                       __global float *b,                       \n" \
"                       __global float *c,                       \n" \
"                       const int n)                          \n" \
"{                                                               \n" \
"    int id = get_global_id(0);                                  \n" \
"    if (id < n)                                                 \n" \
"        c[id] = a[id] + b[id];                                  \n" \
"}                                                               \n";

double vecAddCLMemOnce(cl_kernel *vecAdd,
                       cl_mem *cl_a,
                       cl_mem *cl_b,
                       cl_mem *cl_c,
                       size_t dataLength,
                       cl_command_queue *exeCmdQueue) {

    cl_wrapper::checkError(clSetKernelArg(*vecAdd, 0, sizeof(*cl_a), cl_a));
    cl_wrapper::checkError(clSetKernelArg(*vecAdd, 1, sizeof(*cl_b), cl_b));
    cl_wrapper::checkError(clSetKernelArg(*vecAdd, 2, sizeof(*cl_c), cl_c));
    const auto cl_dataLength = static_cast<cl_int> (dataLength);
    cl_wrapper::checkError(clSetKernelArg(*vecAdd, 3, sizeof(cl_dataLength), &cl_dataLength));

    const size_t global_size[] = {dataLength};
    const size_t wg_size[] = {1024};
    cl_event e;
    cl_wrapper::checkError(clEnqueueNDRangeKernel(*exeCmdQueue,
                                                  *vecAdd,
                                                  1, nullptr,
                                                  global_size,
                                                  wg_size, 0, nullptr, &e));
    clWaitForEvents(1, &e);
    return cl_wrapper::getExecutionTime(&e);
}

double vecAddZeroCopy(cl_wrapper *cl,
                      cl_kernel *vecAdd,
                      ZeroCopyMem<float> *a,
                      ZeroCopyMem<float> *b,
                      ZeroCopyMem<float> *c,
                      size_t dataLength,
                      cl_command_queue *exeCmdQueue) {

    //Measure time include the time for applying zero copy
    ZeroCopyMem<float> tmp_a = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> tmp_b = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> tmp_c = init_zero_copy_region<float>(cl, dataLength);

    double executionTime = vecAddCLMemOnce(vecAdd,
                                           &a->deviceBuffer,
                                           &b->deviceBuffer,
                                           &c->deviceBuffer,
                                           dataLength,
                                           exeCmdQueue);

    clReleaseMemObject(tmp_a.deviceBuffer);
    clReleaseMemObject(tmp_b.deviceBuffer);
    clReleaseMemObject(tmp_c.deviceBuffer);
    return executionTime;
}


double vecAddNonZeroCopy(cl_wrapper *cl,
                         cl_kernel *vecAdd,
                         const float *cpu_a,
                         const float *cpu_b,
                         float *gpu_c,
                         size_t dataLength,
                         cl_command_queue *exeCmdQueue) {
    cl_mem cl_a = cl_make_array(cl, cpu_a, dataLength);
    cl_mem cl_b = cl_make_array(cl, cpu_b, dataLength);
    cl_mem cl_c = cl_make_array<float>(cl, 0, dataLength);

    double executionTime = vecAddCLMemOnce(vecAdd,
                                           &cl_a, &cl_b, &cl_c,
                                           dataLength,
                                           exeCmdQueue);

    cl_pull_array(cl, &cl_c, gpu_c, dataLength);
    clReleaseMemObject(cl_a);
    clReleaseMemObject(cl_b);
    clReleaseMemObject(cl_c);
    return executionTime;
}

int main(int argc, char **argv) {
    //Experiment config
    const size_t LOOP_TIME = 100;
    const size_t w = 1024;
    const size_t h = 1024;
    const size_t dataLength = w * h * 512 / sizeof(float);
    //-Experiment config

    //OpenCL init
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    auto *cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices->deviceId);
    cl_program program = cl->createProgram(&kernelSource, 1);
    cl_kernel vecAdd = cl->createKernel("vecAdd", program);
    cl_command_queue exeCmdQueue = cl->createProfilingCmdQueue();
    //-OpenCL init

    ZeroCopyMem<float> a = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> b = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> c = init_zero_copy_region<float>(cl, dataLength);
    auto *p_a = a.hostPtr;
    auto *p_b = b.hostPtr;
    auto *gpu_c = c.hostPtr;
    auto *cpu_c = (float *) malloc(sizeof(float) * dataLength);
    randArrayGenerator<float>(0.0f, 1.0, p_a, dataLength);
    randArrayGenerator<float>(0.0f, 1.0, p_b, dataLength);
    //memset(gpu_c, 0, dataLength * sizeof(float));
    //memset(cpu_c, 0, dataLength * sizeof(float));
    //warm-up run
    vecAddCLMemOnce(&vecAdd, &a.deviceBuffer, &b.deviceBuffer, &c.deviceBuffer, dataLength, &exeCmdQueue);
    cl_event e;
    cl_int cl_zero = 0;
    cl_wrapper::checkError(
            clEnqueueFillBuffer(exeCmdQueue, c.deviceBuffer, &cl_zero, sizeof(cl_int),
                                0, dataLength * sizeof(cl_int), 0, nullptr, &e));
    clWaitForEvents(1, &e);

    cl_wrapper::checkError(clSetKernelArg(vecAdd, 0, sizeof(a.deviceBuffer), &a.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(vecAdd, 1, sizeof(b.deviceBuffer), &b.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(vecAdd, 2, sizeof(c.deviceBuffer), &c.deviceBuffer));
    const auto cl_dataLength = static_cast<cl_int> (dataLength);
    cl_wrapper::checkError(clSetKernelArg(vecAdd, 3, sizeof(cl_dataLength), &cl_dataLength));

    double t = getCurrentTime();
    for (int j = 0; j < LOOP_TIME; ++j) {
        vecAddCPU(a.hostPtr, b.hostPtr, cpu_c, dataLength);
    }
    flushed_printf("Data Length: %d, CPU exe time %f sec \n",
                   dataLength, (getCurrentTime() - t) / LOOP_TIME);

    int numParts;
    double profTime = 0;
    int numCLQ = 4;
    cl_command_queue exeCmdQueues[numCLQ];
    cl_event exeEvents[numCLQ];
    for (int l = 0; l < numCLQ; ++l) {
        exeCmdQueues[l] = cl->createProfilingCmdQueue();
    }
//    for (int k = 0; k < 13; ++k) { //real benchmark
    for (int k = 0; k < 1; ++k) {
        numParts = pow2(k);
        flushed_printf("-- Running, numParts=%d ...\n", numParts);
        size_t workLength = dataLength / numParts;
        size_t global_size[] = {workLength};
        size_t wg_size[] = {1024};

        for (int m = 1; m <= numCLQ; ++m) {
            flushed_printf("NumQ: %d, ", m);
            t = getCurrentTime();
            for (int j = 0; j < LOOP_TIME; ++j) {
                for (int i = 0; i < numParts; i += m) {
                    for (int l = 0; l < m && l + i < numParts; ++l) {
                        size_t offset[] = {static_cast<size_t>(workLength * (i + l))};
                        cl_wrapper::checkError(clEnqueueNDRangeKernel(
                                exeCmdQueues[l],
                                vecAdd,
                                1, offset,
                                global_size,
                                wg_size, 0, nullptr, &exeEvents[l]));
                    }
                    for (int l = 0; l < m && l + i < numParts; ++l) {
                        clWaitForEvents(1, &exeEvents[l]);
                        profTime += cl_wrapper::getExecutionTime(&exeEvents[l]);
                    }
//                size_t offset[] = {static_cast<size_t>(workLength * i)};
//                cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueue,
//                                                              vecAdd,
//                                                              1, offset,
//                                                              global_size,
//                                                              wg_size, 0, nullptr, &e));
//                clWaitForEvents(1, &e);
//                profTime += cl_wrapper::getExecutionTime(&e);
                }
            }
            flushed_printf("Exe time %f sec, prof run time %f ",
                           (getCurrentTime() - t) / LOOP_TIME, profTime / LOOP_TIME);
            flushed_printf("GPU_CPU diff: %f \n", calVecDiff(gpu_c, cpu_c, dataLength));
            cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c.deviceBuffer, &cl_zero, sizeof(cl_int),
                                                       0, dataLength * sizeof(cl_int), 0, nullptr, &e));
            clWaitForEvents(1, &e);
            profTime = 0;
        }
    }

    printMatrix(p_a, w, h);
    printMatrix(p_b, w, h);
    t = getCurrentTime();
    for (int i = 0; i < LOOP_TIME; ++i) {
        vecAddCPU(p_a, p_b, cpu_c, dataLength);
    }
    flushed_printf("Compute time: %lf sec \n\n", getCurrentTime() - t);
//
//    Zero copy, apply for zero copy regions every time
    memset(c.hostPtr, 0, dataLength * sizeof(c.hostPtr[0]));
    vecAddZeroCopy(cl, &vecAdd, &a, &b, &c, dataLength, &exeCmdQueue);
    printf("Verification zero copy: %f\n", calVecDiff(cpu_c, c.hostPtr, dataLength));
//    for (int i = 0; i < dataLength; ++i) {
//        flushed_printf("c[%d] = %f\n", i, c.hostPtr[i]);
//    }
    memset(gpu_c, 0, dataLength * sizeof(*gpu_c));
    vecAddNonZeroCopy(cl, &vecAdd, p_a, p_b, gpu_c, dataLength, &exeCmdQueue);
    printf("Verification non zero copy: %f\n\n", calVecDiff(cpu_c, gpu_c, dataLength));

    printf("Zero copy, apply for zero copy regions...\n");
    double timer = getCurrentTime();
    for (int i = 0; i < LOOP_TIME; ++i) {
        vecAddZeroCopy(cl, &vecAdd, &a, &b, &c, dataLength, &exeCmdQueue);
    }
    printf("time: %f sec \n\n\n", getCurrentTime() - timer);
    //-Zero copy

    //Zero copy, pure
    printf("Zero copy, pure...\n");
    timer = getCurrentTime();
    for (int i = 0; i < LOOP_TIME; ++i) {
        vecAddCLMemOnce(&vecAdd, &a.deviceBuffer, &b.deviceBuffer, &c.deviceBuffer, dataLength, &exeCmdQueue);
    }
    printf("Zero copy time: %f sec \n\n\n", getCurrentTime() - timer);
    //-Zero copy

    // non zero-copy, copy on every kernel call
    printf("Non zero-copy, copy on every kernel call...\n");
    timer = getCurrentTime();
    for (int i = 0; i < LOOP_TIME; ++i) {
        vecAddNonZeroCopy(cl, &vecAdd, p_a, p_b, gpu_c, dataLength, &exeCmdQueue);
    }
    printf("Non zero copy time: %f sec \n\n\n", getCurrentTime() - timer);
    //-non zero-copy, copy on every kernel call

    //non zero copy, copy once
    printf("Non zero-copy, copy once...\n");
    cl_mem cl_a = cl_make_array(cl, p_a, dataLength);
    cl_mem cl_b = cl_make_array(cl, p_b, dataLength);
    cl_mem cl_c = cl_make_array<float>(cl, 0, dataLength);
    timer = getCurrentTime();
    for (int i = 0; i < LOOP_TIME; ++i) {
        vecAddCLMemOnce(&vecAdd, &cl_a, &cl_b, &cl_c, dataLength, &exeCmdQueue);
    }
    printf("Non zero copy time: %f sec \n\n\n", getCurrentTime() - timer);
    // non zero copy, copy once

    clReleaseMemObject(cl_a);
    clReleaseMemObject(cl_b);
    clReleaseMemObject(cl_c);
    clReleaseEvent(e);
    clReleaseMemObject(a.deviceBuffer);
    clReleaseMemObject(b.deviceBuffer);
    clReleaseMemObject(c.deviceBuffer);
    free(cl);


    return 0;
}
