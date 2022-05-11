//
// Created by pfxu on 8/6/19.
//

#include "heteroCompLib.h"
#include <cstdio>

//const char* matMulCLSrc[] = {
const char* KERNEL_SOURCE[] = {
        "__kernel void mmul(                            \n"
        "const int N,                                   \n"
        "__global float* A,                             \n"
        "__global float* B,                             \n"
        "__global float* C)                             \n"
        "{                                              \n"
        "    int k, j;                                  \n"
        "    int i = get_global_id(0);                  \n"
        "    float tmp;                                 \n"
        "    if (i < N) {                               \n"
        "        for (j = 0; j < N; j++) {              \n"
        "            tmp = 0.0f;                        \n"
        "            for (k = 0; k < N; k++)            \n"
        "                tmp += A[i*N+k] * B[k*N+j];    \n"
        "            C[i*N+j] = tmp;                    \n"
        "        }                                      \n"
        "    }                                          \n"
        "}                                              \n"
};

using namespace std;

static const cl_uint KERNEL_SOURCE_LEN = sizeof(KERNEL_SOURCE) / sizeof(const char *);

int main(int argc, char **argv) {
    //int N = 64;
    //int N = 64 * 64;
    int matDim = 64;
    int N = matDim * matDim;
    std::vector<float> a, b, c;
    a.resize(N);
    b.resize(N);
    c.resize(N);
    cout << "Generating random vectors[" << N << "]... " << endl;
    randArrayGenerator<float>(-10.0, 10.0, a.data(), N);
    randArrayGenerator<float>(-10.0, 10.0, b.data(), N);
    memset(&c[0], 0, N * sizeof(c[0]));

    cl_int error;
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    auto *cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    cl_program program = cl->createProgram(KERNEL_SOURCE, KERNEL_SOURCE_LEN);
    cl_kernel kernel_vecAdd = cl->createKernel("mmul", program);

    cl_mem cl_a = cl_make_array(cl, a.data(), N);
    cl_mem cl_b = cl_make_array(cl, b.data(), N);
    cl_mem cl_c = cl_make_array(cl, c.data(), N);

    cl_wrapper::checkError(clSetKernelArg(kernel_vecAdd, 0, sizeof(N), &N));
    cl_wrapper::checkError(clSetKernelArg(kernel_vecAdd, 1, sizeof(cl_a), &cl_a));
    cl_wrapper::checkError(clSetKernelArg(kernel_vecAdd, 2, sizeof(cl_b), &cl_b));
    cl_wrapper::checkError(clSetKernelArg(kernel_vecAdd, 3, sizeof(cl_c), &cl_c));

    size_t resolution;
    clGetDeviceInfo(cl->device, CL_DEVICE_PROFILING_TIMER_RESOLUTION,
                    sizeof(resolution), &resolution, nullptr);
    cout << "Gathering info, timer resolution: " << resolution << endl;
    cl_uint num_queues = 32;
    cl_command_queue exeCmdQueues[num_queues];
    cl_event kernel_events[num_queues];
    for (auto &exeCmdQueue : exeCmdQueues) exeCmdQueue = cl->createProfilingCmdQueue();
    const size_t local_work_size[2] = {8, 8};
    const size_t global_work_size[2] = {static_cast<size_t>(matDim), static_cast<size_t>(matDim)};

    for (int m = 1; m <= num_queues; m += 2) {
        printf("Number of queues: %d\n", m);
        cl_event user_event = clCreateUserEvent(cl->context, &error);
        cl_wrapper::checkError(error);
        for (int i = 0; i < m; ++i) {
            clEnqueueNDRangeKernel(exeCmdQueues[i], kernel_vecAdd,
                                   1,
                                   nullptr,
                                   global_work_size,
                                   local_work_size,
                                   1, &user_event, &kernel_events[i]);
//                                   0, nullptr, &kernel_events[i]);
        }
        clSetUserEventStatus(user_event, CL_COMPLETE);
        clWaitForEvents(m, kernel_events);
        clReleaseEvent(user_event);

        cl_ulong buff = 0;
        cl_ulong first = 0;
        printf("%-13s%-13s%-13s%-13s%-13s\n", "Queued", "Submit", "Start", "End", "Complete");
        for (int i = 0; i < m; ++i) {
            clGetEventProfilingInfo(kernel_events[i], CL_PROFILING_COMMAND_QUEUED,
                                    sizeof(buff), &buff, nullptr);
            if (i == 0) first = buff;
            printf("%-12lu ", buff - first);
            clGetEventProfilingInfo(kernel_events[i], CL_PROFILING_COMMAND_SUBMIT,
                                    sizeof(buff), &buff, nullptr);
            printf("%-12lu ", buff - first);
            clGetEventProfilingInfo(kernel_events[i], CL_PROFILING_COMMAND_START,
                                    sizeof(buff), &buff, nullptr);
            printf("%-12lu ", buff - first);
            clGetEventProfilingInfo(kernel_events[i], CL_PROFILING_COMMAND_END,
                                    sizeof(buff), &buff, nullptr);
            printf("%-12lu ", buff - first);
            clGetEventProfilingInfo(kernel_events[i], CL_PROFILING_COMMAND_COMPLETE,
                                    sizeof(buff), &buff, nullptr);
            printf("%-12lu\n", buff - first);
        }
    }


    clReleaseMemObject(cl_a);
    clReleaseMemObject(cl_b);
    clReleaseMemObject(cl_c);

    return 0;
}
