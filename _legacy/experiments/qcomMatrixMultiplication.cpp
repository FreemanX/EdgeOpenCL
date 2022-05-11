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
        "        b = vload4(0, matrix_b + j * matrix_b_width + (wid_x * 4));\n", // vload4 load 4 floats
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

struct matrix_t {
    size_t width, height, dataLength;
    float *elements;
};

matrix_t genMatrix(size_t height, size_t width) {
    matrix_t mat{};
    mat.width = width;
    mat.height = height;
    mat.dataLength = width * height;
    mat.elements = (float *) malloc(sizeof(float) * mat.dataLength);
    srand(static_cast<unsigned int>(getCurrentTime()));
    randArrayGenerator(-10.0f, 10.0f, mat.elements, mat.width * mat.height);

    return mat;
}

double gpu_gemm_qcom(const matrix_t *mat_A, const matrix_t *mat_B, matrix_t *mat_C) {
    const size_t mat_A_size = mat_A->dataLength;
    const size_t mat_B_size = mat_B->dataLength;
    const size_t mat_C_size = mat_C->dataLength;

    //Init OpenCL
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    auto *cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    cl_program program = cl->createProgram(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN);
    cl_kernel kernel_8x4 = cl->createKernel("matmul_8x4_blocks", program);
    cl_kernel kernel_rem = cl->createKernel("matmul_remainder", program);
    cl_command_queue exeCmdQueue = cl->createProfilingCmdQueue();

    double curTime = getCurrentTime();

    //Execute 1st kernel kernel_8x4
    cl_mem mat_a_mem = cl_make_array<float const>(cl, mat_A->elements, mat_A_size);
    cl_mem mat_b_mem = cl_make_array<float const>(cl, mat_B->elements, mat_B_size);
    cl_mem mat_c_mem = cl_make_array<float>(cl, mat_C->elements, mat_C_size); 
    const auto matrix_b_width = static_cast<const cl_int>(mat_B->width);
    const auto matrix_a_width = static_cast<const cl_int>(mat_A->width);
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 0, sizeof(mat_a_mem), &mat_a_mem));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 1, sizeof(mat_b_mem), &mat_b_mem));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 2, sizeof(mat_c_mem), &mat_c_mem));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 3, sizeof(matrix_b_width), &matrix_b_width));
    cl_wrapper::checkError(clSetKernelArg(kernel_8x4, 4, sizeof(matrix_a_width), &matrix_a_width));

    const size_t tiled_global_work_size[] = {static_cast<size_t>(mat_B->width / 4),
                                             static_cast<size_t>(mat_A->height / 8)};
    if (tiled_global_work_size[0] != 0 && tiled_global_work_size[1] != 0) {
        cl_event e;
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
    const auto x_rem_start = static_cast<const cl_int>((mat_B->width / 4) * 4);
    const auto y_rem_start = static_cast<const cl_int>((mat_A->height / 8) * 8);
    const cl_int right_y_rem_start = 0;
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 0, sizeof(mat_a_mem), &mat_a_mem));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 1, sizeof(mat_b_mem), &mat_b_mem));
    cl_wrapper::checkError(clSetKernelArg(kernel_rem, 2, sizeof(mat_c_mem), &mat_c_mem));
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
                nullptr, 0, nullptr, nullptr));
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
                nullptr, 0, nullptr, nullptr);
    }

    //Get result
    cl_pull_array(cl, &mat_c_mem, mat_C->elements, mat_C_size);

    curTime = getCurrentTime() - curTime;

    // Clean up cl resources that aren't automatically handled by cl_wrapper
    clReleaseMemObject(mat_a_mem);
    clReleaseMemObject(mat_b_mem);
    clReleaseMemObject(mat_c_mem);

    return curTime;
}

int main(int argc, char **argv) {
    //const size_t m = 32;
    //const size_t n = 43264;
    //const size_t k = 144;
    size_t m, n, k;
    cout << endl << "m: ";
    cin >> m;
    cout << endl << "n: ";
    cin >> n;
    cout << endl << "k: ";
    cin >> k;

    //Init matrix
    matrix_t mat_A = genMatrix(m, k);
    matrix_t mat_B = genMatrix(k, n);
    matrix_t mat_C_gpu{};
    mat_C_gpu.height = mat_A.height;
    mat_C_gpu.width = mat_B.width;
    mat_C_gpu.dataLength = mat_C_gpu.width * mat_C_gpu.height;
    mat_C_gpu.elements = (float*) malloc(mat_C_gpu.height * mat_C_gpu.width * sizeof(float));
    memset(&mat_C_gpu.elements[0], 0, mat_C_gpu.dataLength * sizeof(mat_C_gpu.elements[0]));

    printf("A(%ld, %ld) x B(%ld, %ld) = C(%ld, %ld)\n",
           mat_A.height, mat_A.width, mat_B.height, mat_B.width, mat_C_gpu.height, mat_C_gpu.width);
    std::cout << "A: " << std::endl;
    printMatrix(mat_A.elements, static_cast<int>(mat_A.height), static_cast<int>(mat_A.width));
    std::cout << "B: " << std::endl;
    printMatrix(mat_B.elements, static_cast<int>(mat_B.height), static_cast<int>(mat_B.width));

    int repeat = 1;
    double runTimeGPU = 0.0;
    for (int i = 0; i < repeat; ++i) {
        runTimeGPU += gpu_gemm_qcom(&mat_A, &mat_B, &mat_C_gpu);
    }

    cout << "GPU Run time: " << runTimeGPU / repeat << "s" << endl;

    // CPU
    matrix_t mat_C_cpu{};
    mat_C_cpu.height = mat_A.height;
    mat_C_cpu.width = mat_B.width;
    mat_C_cpu.dataLength = mat_C_cpu.height * mat_C_cpu.width;
    mat_C_cpu.elements = (float *) malloc(sizeof(float) * mat_C_cpu.dataLength);
    memset(&mat_C_cpu.elements[0], 0, mat_C_cpu.dataLength * sizeof(mat_C_cpu.elements[0]));

    double runtimeCPU = 0.0;
    for (int i = 0; i < repeat; ++i) {
        memset(&mat_C_cpu.elements[0], 0, mat_C_cpu.dataLength * sizeof(mat_C_cpu.elements[0]));
        runtimeCPU += gemm_cpu(0, 0, m, n, k, 1,
                               mat_A.elements, static_cast<int>(mat_A.width),
                               mat_B.elements, static_cast<int>(mat_B.width),
                               1,
                               mat_C_cpu.elements, static_cast<int>(mat_C_cpu.width));
    }

    std::cout << "CPU Run time: " << runtimeCPU / repeat << "s" << std::endl;
    std::cout << "CPU/GPU: " << runtimeCPU / runTimeGPU << std::endl;

    std::cout << "Diff: " << calVecDiff(mat_C_cpu.elements, mat_C_gpu.elements, m * n) << "\n";
    std::cout << "CPU:" << std::endl;
    printMatrix(mat_C_cpu.elements, static_cast<int>(mat_C_cpu.height), static_cast<int>(mat_C_cpu.width));
    std::cout << "GPU:" << std::endl;
    printMatrix(mat_C_gpu.elements, static_cast<int>(mat_C_cpu.height), static_cast<int>(mat_C_cpu.width));

    return 0;
}















