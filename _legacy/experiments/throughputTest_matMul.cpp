//
// Created by pfxu on 8/13/19.
//

#include "heteroCompLib.h"
#include "CPU_blas.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>

using namespace std;

static const char *PROGRAM_SOURCE[] = {
        "__kernel void matmul_8x4_blocks(__global const float *matrix_a,\n",
        "                                __global const float *matrix_b,\n",
        "                                __global       float *matrix_c,\n",
        "                                               int    matrix_b_width,\n",
        "                                               int    matrix_a_width)\n",
        "{\n",
        "    const int wid_x = get_global_id(0);\n",
        "    const int wid_y = get_global_id(1);\n",
        "\n",
        "    float  a[8];\n",
        "    float4 b;\n",
        "    float4 c[8];\n",
        "\n",
        "    for (int i = 0; i < 8; ++i)\n",
        "    {\n",
        "        c[i] = (float4)(0.0f);\n",
        "    }\n",
        "\n",
        "    for (int j = 0; j < matrix_a_width; ++j)\n",
        "    {\n",
        "        b = vload4(0, matrix_b + j * matrix_b_width + (wid_x * 4));\n",
        //"        b = 1.0;\n",
        "\n",
        "#pragma unroll\n",
        "        for (int i = 0; i < 8; ++i)\n",
        "        {\n",
        "            a[i] = matrix_a[((wid_y * 8) + i) * matrix_a_width + j];\n",
        "        }\n",
        "\n",
        "#pragma unroll\n",
        "        for (int i = 0; i < 8; ++i)\n",
        "        {\n",
        "            c[i] += a[i] * b;\n",
        "        }\n",
        "    }\n",
        "\n",
        "#pragma unroll\n",
        "    for (int i = 0; i < 8; ++i)\n",
        "    {\n",
        "        vstore4(c[i], 0, matrix_c + ((wid_y * 8) + i) * matrix_b_width + (wid_x * 4));\n",
        "    }\n",
        "}\n",
        "\n",
        "__kernel void matmul_remainder(__global const  float *matrix_a,\n",
        "                               __global const  float *matrix_b,\n",
        "                               __global        float *matrix_c,\n",
        "                                               int    x_rem_start,\n",
        "                                               int    y_rem_start,\n",
        "                                               int    matrix_b_width,\n",
        "                                               int    matrix_a_width)\n",
        "{\n",
        "    const int wid_x = get_global_id(0) + x_rem_start;\n",
        "    const int wid_y = get_global_id(1) + y_rem_start;\n",
        "\n",
        "    float c     = 0.0f;\n",
        "    int   a_idx = matrix_a_width * wid_y;\n",
        "    int   b_idx = wid_x;\n",
        "\n",
        "#pragma unroll 8\n",
        "    for (int i = 0; i < matrix_a_width; ++i)\n",
        "    {\n",
        "        c += matrix_a[a_idx] * matrix_b[b_idx];\n",
        "        ++a_idx;\n",
        "        b_idx += matrix_b_width;\n",
        "    }\n",
        "\n",
        "    const int c_idx = wid_x + matrix_b_width * wid_y;\n",
        "    matrix_c[c_idx] = c;\n",
        "}\n"
};

static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_kernel kernel_8x4;
cl_kernel kernel_rem;
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

double gpu_gemm(matrix_t *mat_A,
                matrix_t *mat_B,
                matrix_t *mat_C) {

    double runtimeGPU = getCurrentTime();
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 0, sizeof(cl_mem), &mat_A->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 1, sizeof(cl_mem), &mat_B->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 2, sizeof(cl_mem), &mat_C->data.deviceBuffer));

    const auto matrix_b_width = static_cast<const cl_int>(mat_B->width);
    const auto matrix_a_width = static_cast<const cl_int>(mat_A->width);
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 3, sizeof(matrix_b_width), &matrix_b_width));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 4, sizeof(matrix_a_width), &matrix_a_width));

    const size_t tiled_global_work_size[] = {static_cast<size_t>(mat_B->width / 4),
                                             static_cast<size_t>(mat_A->height / 8)};
    cl_event e;
    if (tiled_global_work_size[0] != 0 && tiled_global_work_size[1] != 0) {
        cl_wrapper::checkError(
                clEnqueueNDRangeKernel(
                        exeCmdQueue,
                        kernel_8x4,
                        2, nullptr,
                        tiled_global_work_size,
                        nullptr, 0, nullptr, &e));
        clWaitForEvents(1, &e);
    }

    //Execute 2nd kernel
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 0, sizeof(cl_mem), &mat_A->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 1, sizeof(cl_mem), &mat_B->data.deviceBuffer));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 2, sizeof(cl_mem), &mat_C->data.deviceBuffer));

    const auto x_rem_start = static_cast<const cl_int>((mat_B->width / 4) * 4);
    const auto y_rem_start = static_cast<const cl_int>((mat_A->height / 8) * 8);
    const cl_int right_y_rem_start = 0;
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 3, sizeof(x_rem_start), &x_rem_start));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 4, sizeof(right_y_rem_start), &right_y_rem_start));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 5, sizeof(matrix_b_width), &matrix_b_width));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 6, sizeof(matrix_a_width), &matrix_a_width));

    const size_t right_rem_work_size[] = {static_cast<size_t>(mat_B->width - x_rem_start),
                                          static_cast<size_t>(mat_A->height)};
    if (right_rem_work_size[0] != 0 && right_rem_work_size[1] != 0) {
        cl_wrapper::checkError(clEnqueueNDRangeKernel(
                exeCmdQueue,
                kernel_rem,
                2, nullptr,
                right_rem_work_size,
                nullptr, 0, nullptr, &e));
        clWaitForEvents(1, &e);
    }

    const cl_int bottom_x_rem_start = 0;
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 3, sizeof(bottom_x_rem_start), &bottom_x_rem_start));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 4, sizeof(y_rem_start), &y_rem_start));
    //Covers the remaining bottom portion of the result matrix not covered above.
    const size_t bottom_rem_work_size[] = {static_cast<size_t>(x_rem_start),
                                           static_cast<size_t>(mat_A->height - y_rem_start)};
    if (bottom_rem_work_size[0] != 0 && bottom_rem_work_size[1] != 0) {
        clEnqueueNDRangeKernel(
                exeCmdQueue,
                kernel_rem,
                2, nullptr,
                bottom_rem_work_size,
                nullptr, 0, nullptr, &e);
        clWaitForEvents(1, &e);
    }

    //Get result
    runtimeGPU = getCurrentTime() - runtimeGPU;
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
    size_t LOOP_TIME = 10000;

    //Init CL
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    program = cl->createProgram(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN);
    kernel_8x4 = cl->createKernel("matmul_8x4_blocks", program);
    kernel_rem = cl->createKernel("matmul_remainder", program);
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


        setCurThreadAffinity(bigCPUList);
        setCurThreadAffinity(bigCPUList);
        double runtimeGPU = 0;
        for (int i = 0; i < LOOP_TIME; ++i) {
            runtimeGPU += gpu_gemm(&mat_A, &mat_B, &mat_C);
        }
        std::cout << "Size = " << size << " GPU Run time: " << runtimeGPU / LOOP_TIME << " s" << std::endl;

        setCurThreadAffinity(7);
        setCurThreadAffinity(7);
        omp_set_num_threads(1);
        double runtimeCPU = 0.0;
        for (int i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += gemm_cpu(
                    0, 0, m, n, k, 1,
                    mat_A.data.hostPtr, static_cast<int>(mat_A.width),
                    mat_B.data.hostPtr, static_cast<int>(mat_B.width),
                    1,
                    cpu_c, static_cast<int>(mat_C.width));
        }
        std::cout << "Size = " << size  << " Big1 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;

        setCurThreadAffinity(3);
        setCurThreadAffinity(3);
        omp_set_num_threads(1);
        runtimeCPU = 0.0;
        for (int i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += gemm_cpu(
                    0, 0, m, n, k, 1,
                    mat_A.data.hostPtr, static_cast<int>(mat_A.width),
                    mat_B.data.hostPtr, static_cast<int>(mat_B.width),
                    1,
                    cpu_c, static_cast<int>(mat_C.width));
        }
        std::cout << "Size = " << size  << " Little1 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;

        setCurThreadAffinity(bigCPUList);
        setCurThreadAffinity(bigCPUList);
        omp_set_num_threads(4);
        runtimeCPU = 0.0;
        for (int i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += gemm_cpu(
                    0, 0, m, n, k, 1,
                    mat_A.data.hostPtr, static_cast<int>(mat_A.width),
                    mat_B.data.hostPtr, static_cast<int>(mat_B.width),
                    1,
                    cpu_c, static_cast<int>(mat_C.width));
        }
        std::cout << "Size = " << size  << " Big4 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;

        setCurThreadAffinity(littleCPUList);
        setCurThreadAffinity(littleCPUList);
        omp_set_num_threads(4);
        runtimeCPU = 0.0;
        for (int i = 0; i < LOOP_TIME; ++i) {
            memset(cpu_c, 0, mat_C.dataLength * sizeof(float));
            runtimeCPU += gemm_cpu(
                    0, 0, m, n, k, 1,
                    mat_A.data.hostPtr, static_cast<int>(mat_A.width),
                    mat_B.data.hostPtr, static_cast<int>(mat_B.width),
                    1,
                    cpu_c, static_cast<int>(mat_C.width));
        }
        std::cout << "Size = " << size  << " Little4 Run time: " << runtimeCPU / LOOP_TIME << " s" << std::endl;


//        cout << "Diff: " << calVecDiff(cpu_c, mat_C.data.hostPtr, mat_C.dataLength) << endl;

        // Clean up cl resources that aren't automatically handled by cl_wrapper
        clReleaseMemObject(mat_A.data.deviceBuffer);
        clReleaseMemObject(mat_A.data.deviceBuffer);
        clReleaseMemObject(mat_A.data.deviceBuffer);
        LOOP_TIME = LOOP_TIME / 10 > 10 ? LOOP_TIME / 10 : 10;
    }
    return 0;
}
