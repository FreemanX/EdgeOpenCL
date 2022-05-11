//
// Created by pfxu on 3/2/19.
//
#include "heteroCompLib.h"
#include <cstdio>


static const char *KERNEL_SOURCE[] = {
        "__kernel void vecAdd(                          \n"
        "const int N,                                   \n"
        "__global float* A,                             \n"
        "__global float* B,                             \n"
        "__global float* C)                             \n"
        "{                                              \n"
        "    int i = get_global_id(0);                  \n"
        "    if (i < N) {                               \n"
        "       C[i] = A[i] + B[i];                     \n"
//                        "       float tmp[20];                                 \n"
//                        "       for(int j = 0; j < 10000; j++){              \n" // meaningless loop for extend exe time
//                        "           for(int n = 0; n < 20; n++) {                   \n"
//                        "               tmp[n] = (B[j%N] / A[j%N])/1000000.0;   \n"
//                        "           }                                          \n"
//                        "       }                                              \n"
//                        "       for(int j = 0; j < 20; j++) {           \n"
//                        "           C[i] = C[i] + tmp[j];               \n"
//                        "       }                                       \n"
        "    }                                          \n"
        "}                                              \n"
};

using namespace std;

static const cl_uint KERNEL_SOURCE_LEN = sizeof(KERNEL_SOURCE) / sizeof(const char *);

int main(int argc, char **argv) {
    //int N = 64;
    //int N = 64 * 64;
    int N = 1024 * 1024 * 512 / sizeof(float);
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
    cl_kernel kernel_vecAdd = cl->createKernel("vecAdd", program);

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
//    const size_t local_work_size[] = {1024};

    for (int m = 32; m <= num_queues; m += 2) {
        printf("Number of queues: %d\n", m);
        cl_event user_event = clCreateUserEvent(cl->context, &error);
        cl_wrapper::checkError(error);
        for (int i = 0; i < m; ++i) {
            size_t offset[] = {static_cast<size_t>(i * (N / m))};
            size_t global_work_size[] = {static_cast<size_t>(N / m)};
            clEnqueueNDRangeKernel(exeCmdQueues[i], kernel_vecAdd,
                                   1,
                                   offset,
                                   global_work_size,
                                   nullptr,
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
