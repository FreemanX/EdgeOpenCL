//
// Created by pfxu on 21/8/2019.
//

#include "heteroCompLib.h"
#include "CPU_blas.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>

using namespace std;

static const char *PROGRAM_SOURCE[] = {
        "__kernel void matAdd(__global const float* A,\n",
        "	__global const float* B,\n",
        "	__global float* C,\n",
        "	const int N) {\n",
        "	int id = get_global_id(0);\n",
        "	int idx = mul24(id, 4);\n",
        "	float4 tmpA = vload4(idx, A);\n",
        "	float4 tmpB = vload4(idx, B);\n",
        "	(*((__global float4*)&C[idx])) = (tmpA+tmpB);\n",
        "}\n",
        "\n",
        "__kernel void matAdd2(\n",
        "	__global float* matrixA,\n",
        "	__global float* matrixB,\n",
        "	__global float* MatrixSum,\n",
        "	const int rows,\n",
        "	const int cols)\n",
        "{\n",
        "	int i = get_global_id(0);\n",
        "	int j = get_global_id(1);\n",
        "	int offset = mul24(j, cols);\n",
        "	int index = mad24(i, 4, offset);\n",
        "\n",
        "	float4 tmpA = (*((__global float4*)&matrixA[index]));\n",
        "	float4 tmpB = (*((__global float4*)&matrixB[index]));\n",
        "	(*((__global float4*)&MatrixSum[index])) = (tmpA+tmpB);\n",
        "\n",
        "}\n"
};

static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_kernel kernel_add;
cl_kernel kernel_add2;
cl_command_queue exeCmdQueue;

struct matrix_t {
    size_t width, height, dataLength;
    ZeroCopyMem<float> data;

//    float *data;
    void setSize(size_t size) {
        width = size;
        height = size;
        dataLength = width * height;
        data.size = dataLength * sizeof(float);
        data.length = dataLength;
    }
};

double gpu_matAdd(matrix_t *mat_A, matrix_t *mat_B, matrix_t *mat_C) {
    double runtimeGPU = getCurrentTime();
    const auto N = static_cast<const cl_int>(mat_B->width) * static_cast<const cl_int>(mat_A->width);
    cl_wrapper::checkError(clSetKernelArg(kernel_add, 0, sizeof(cl_mem), &mat_A->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add, 1, sizeof(cl_mem), &mat_B->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add, 2, sizeof(cl_mem), &mat_C->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add, 3, sizeof(N), &N));

    const size_t global_work_size[] = {static_cast<size_t>(N / 4)};
    const size_t local_work_size[] = {1024};
    cl_event event;
    cl_wrapper::checkError(clEnqueueNDRangeKernel(
            exeCmdQueue,
            kernel_add,
            1, nullptr,
            global_work_size,
            local_work_size, 0, nullptr, &event));
    clWaitForEvents(1, &event);

    //Get result
    runtimeGPU = getCurrentTime() - runtimeGPU;
    return runtimeGPU;
}

double gpu_matAdd2(matrix_t *mat_A, matrix_t *mat_B, matrix_t *mat_C) {
    double runtimeGPU = getCurrentTime();
    cl_int width = mat_C->width;
    cl_int height = mat_C->height;
    cl_wrapper::checkError(clSetKernelArg(kernel_add2, 0, sizeof(cl_mem), &mat_A->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add2, 1, sizeof(cl_mem), &mat_B->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add2, 2, sizeof(cl_mem), &mat_C->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_add2, 3, sizeof(cl_int), &height));
    cl_wrapper::checkError(clSetKernelArg(kernel_add2, 4, sizeof(cl_int), &width));

    const size_t global_work_size[] = {static_cast<size_t>(width / 4), static_cast<size_t>(height)};
    size_t localDim = MIN(1024, mat_C->width);
//    size_t localDim = 512;
    const size_t local_work_size[] = {localDim, localDim};
    cl_event event;
    cl_wrapper::checkError(clEnqueueNDRangeKernel(
            exeCmdQueue,
            kernel_add2,
            1, nullptr,
            global_work_size,
            local_work_size,
//            nullptr,
            0, nullptr, &event));
    clWaitForEvents(1, &event);
    //Get result
    runtimeGPU = getCurrentTime() - runtimeGPU;
//    runtimeGPU = cl_wrapper::getExecutionTime(&event);
    return runtimeGPU;
}

int main(int argc, char **argv) {
    cout << "Initiating...\n";
    std::vector<int> littleCPUList;
    std::vector<int> bigCPUList;
    for (int i = 0; i < 4; ++i) {
        littleCPUList.push_back(i);
        bigCPUList.push_back(i + 4);
    }
    //Experiment config
    size_t LOOP_TIME = 1000000;

    //Init CL
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    program = cl->createProgram(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN);
    kernel_add = cl->createKernel("matAdd", program);
    kernel_add2 = cl->createKernel("matAdd2", program);
    exeCmdQueue = cl->createProfilingCmdQueue();
    //-Init CL

    int size = 2;
    for (int q = 0; q < 12; ++q) {
        size = size << 1;
        cout << size << endl;

        size_t m, n, k;
        m = n = k = size;
        matrix_t mat_A{};
        matrix_t mat_B{};
        matrix_t mat_C{};
        mat_A.height = m;
        mat_B.height = k;
        mat_C.height = mat_A.height;
        mat_A.width = k;
        mat_B.width = n;
        mat_C.width = mat_B.width;
        mat_A.dataLength = m * k;
        mat_B.dataLength = n * k;
        mat_C.dataLength = m * n;
        mat_A.data = init_zero_copy_region<float>(cl, mat_A.dataLength);
        mat_B.data = init_zero_copy_region<float>(cl, mat_B.dataLength);
        mat_C.data = init_zero_copy_region<float>(cl, mat_C.dataLength);
        randArrayGenerator<float>(0.0f, 1.0, mat_A.data.hostPtr, mat_A.dataLength);
        randArrayGenerator<float>(0.0f, 1.0, mat_B.data.hostPtr, mat_B.dataLength);
        auto *cpu_c = (float *) malloc(sizeof(float) * mat_C.dataLength);
        memset(cpu_c, 0, mat_C.dataLength * sizeof(float));


        double runtimeGPU = 0;
        setCurThreadAffinity(bigCPUList);
        setCurThreadAffinity(bigCPUList);
        for (size_t i = 0; i < LOOP_TIME; ++i) {
            runtimeGPU += gpu_matAdd2(&mat_A, &mat_B, &mat_C);
        }
        std::cout << "Size = " << size << " GPU Run time: " << runtimeGPU / LOOP_TIME << " s" << std::endl;
//        std::cout << "Size = " << size << " GPU Run time: " << runtimeGPU << " s" << std::endl;

        double runtimeCPU = 0.0;
        setCurThreadAffinity(7);
        setCurThreadAffinity(7);
        omp_set_num_threads(1);
        for (size_t i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += vecAddCPU(mat_A.data.hostPtr, mat_B.data.hostPtr, cpu_c, mat_C.dataLength);
        }
        std::cout << "Size = " << size << " Big1 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;
//        std::cout << "Size = " << size << " Big1 Run time: " << runtimeCPU << " s" << std::endl;

        setCurThreadAffinity(3);
        setCurThreadAffinity(3);
        omp_set_num_threads(1);
        runtimeCPU = 0.0;
        for (size_t i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += vecAddCPU(mat_A.data.hostPtr, mat_B.data.hostPtr, cpu_c, mat_C.dataLength);
        }
        std::cout << "Size = " << size << " Little1 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;
//        std::cout << "Size = " << size << " Little1 Run time: " << runtimeCPU << " s" << std::endl;

        setCurThreadAffinity(bigCPUList);
        setCurThreadAffinity(bigCPUList);
        omp_set_num_threads(4);
        runtimeCPU = 0.0;
        for (size_t i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += vecAddCPU(mat_A.data.hostPtr, mat_B.data.hostPtr, cpu_c, mat_C.dataLength);
        }
        std::cout << "Size = " << size << " Big4 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;
//        std::cout << "Size = " << size << " Big4 Run time: " << runtimeCPU << " s" << std::endl;

        setCurThreadAffinity(littleCPUList);
        setCurThreadAffinity(littleCPUList);
        omp_set_num_threads(4);
        runtimeCPU = 0.0;
        for (size_t i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += vecAddCPU(mat_A.data.hostPtr, mat_B.data.hostPtr, cpu_c, mat_C.dataLength);
        }
        std::cout << "Size = " << size << " Little4 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;
//        std::cout << "Size = " << size << " Little4 Run time: " << runtimeCPU << " s" << std::endl;


//        cout << "Diff: " << calVecDiff(cpu_c, mat_C.data.hostPtr, MIN(100, size)) << endl;

        // Clean up cl resources that aren't automatically handled by cl_wrapper
        clReleaseMemObject(mat_A.data.deviceBuffer);
        clReleaseMemObject(mat_A.data.deviceBuffer);
        clReleaseMemObject(mat_A.data.deviceBuffer);
        LOOP_TIME = LOOP_TIME / 10 > 100 ? LOOP_TIME / 10 : 100;
    }
    return 0;
}
