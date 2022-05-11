//
// Created by pfxu on 6/27/19.
//

#include "heteroCompLib.h"
#include "HeteroComputeAgent.h"
#include <map>

using namespace std;

const char *kernelSource = "\n"\
"__kernel void writeNum(__global float *a, const float n) { \n"\
"    int id = get_global_id(0); \n"\
"    if (id < n) \n"\
"        a[id] = n;  \n"\
"}  \n";


cl_wrapper *cl;
cl_int err;

void initCLEnv() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices->deviceId);
}

ZeroCopyMem<float> createSubZeroCopyMem(ZeroCopyMem<float> buffer, cl_buffer_region *sub_buffer_region) {
    cout << "Creating zero copy sub-buffer {"
         << sub_buffer_region->origin << ", " << sub_buffer_region->size << "}\n";

    ZeroCopyMem<float> sub_buffer{};
    sub_buffer.length = sub_buffer_region->size / sizeof(float);

    sub_buffer.deviceBuffer = clCreateSubBuffer(
            buffer.deviceBuffer,
            CL_MEM_READ_WRITE,
            CL_BUFFER_CREATE_TYPE_REGION, sub_buffer_region, &err);
    cl_wrapper::checkError(err);

    cout << "Mapping hostPtr\n";
    cl_event event;
    sub_buffer.hostPtr = (float *) clEnqueueMapBuffer(
            cl->memCmdQueue,
            sub_buffer.deviceBuffer, CL_TRUE, CL_MAP_WRITE, 0,
            sub_buffer_region->size,
            0, nullptr, &event, &err);
    clWaitForEvents(1, &event);
    cl_wrapper::checkError(err);

    return sub_buffer;
}

int main(int argc, char **argv) {
    initCLEnv();
    size_t dataLength = 10 * 1024;
    ZeroCopyMem<float> buffer = init_zero_copy_region<float>(cl, dataLength);
    for (int i = 0; i < buffer.length; ++i) {
        buffer.hostPtr[i] = 1;
    }
    printMatrix(buffer.hostPtr, 10, 10);

    cl_buffer_region sub_buffer_region{.origin = 0, .size = 1024};
    ZeroCopyMem<float> sub_buffer = createSubZeroCopyMem(buffer, &sub_buffer_region);
    for (int i = 0; i < sub_buffer.length; ++i) {
        sub_buffer.hostPtr[i] = 2;
    }
    printMatrix(buffer.hostPtr + sub_buffer_region.origin, 10, 10);

    cl_buffer_region sub_buffer_region1{.origin = 1024, .size = 332};
    ZeroCopyMem<float> sub_buffer1 = createSubZeroCopyMem(buffer, &sub_buffer_region1);
    for (int i = 0; i < sub_buffer1.length; ++i) {
        sub_buffer1.hostPtr[i] = 3;
    }
    printMatrix(buffer.hostPtr + sub_buffer_region1.origin / sizeof(float), 10, 10);
    printMatrix(sub_buffer1.hostPtr, 10, 10);

    cl_buffer_region super_buffer_region{.origin = 0, .size = buffer.size};
    ZeroCopyMem<float> super_sub_buffer = createSubZeroCopyMem(buffer, &super_buffer_region);
    for (int i = 0; i < super_sub_buffer.length; ++i) {
        super_sub_buffer.hostPtr[i] = 4;
    }
    printMatrix(buffer.hostPtr + (dataLength - 100), 10, 10);
    printMatrix(super_sub_buffer.hostPtr, 10, 10);

    clReleaseMemObject(sub_buffer.deviceBuffer);
    clReleaseMemObject(sub_buffer1.deviceBuffer);

    cl_buffer_region super_buffer_region1{.origin = 0, .size = buffer.size};
    ZeroCopyMem<float> super_sub_buffer1 = createSubZeroCopyMem(buffer, &super_buffer_region1);
    for (int i = 0; i < super_sub_buffer1.length; ++i) {
        super_sub_buffer1.hostPtr[i] = 5;
    }
    printMatrix(buffer.hostPtr + (dataLength - 100), 10, 10);
    printMatrix(super_sub_buffer1.hostPtr, 10, 10);

    clReleaseMemObject(buffer.deviceBuffer);

    size_t dataSize = BUFFER_SIZE;
    double t = getCurrentTime();
    ZeroCopyMem<float> a = init_zero_copy_region_byte<float>(cl, dataSize);
    ZeroCopyMem<float> b = init_zero_copy_region_byte<float>(cl, dataSize);
    ZeroCopyMem<float> c = init_zero_copy_region_byte<float>(cl, dataSize);
    cout << "buffer init time: " << getCurrentTime() - t << "s \n";


    cl_buffer_region bufferRegion{.origin = 0, .size = dataSize};
    t = getCurrentTime();
    ZeroCopyMem<float> sub_a = createSubZeroCopyMem(a, &bufferRegion);
    ZeroCopyMem<float> sub_b = createSubZeroCopyMem(b, &bufferRegion);
    ZeroCopyMem<float> sub_c = createSubZeroCopyMem(c, &bufferRegion);
    cout << "sub buffer init time: " << getCurrentTime() - t << "s \n";

//    cout << "fill sub buffer\n";
//    randArrayGenerator<float>(0.0, 1.0, sub_a.hostPtr, sub_a.length);

    cout << "Release buffer\n";
    clReleaseMemObject(a.deviceBuffer);
    clReleaseMemObject(b.deviceBuffer);
    clReleaseMemObject(c.deviceBuffer);

    cout << "Release sub-buffer\n";
    clReleaseMemObject(sub_a.deviceBuffer);
    clReleaseMemObject(sub_b.deviceBuffer);
    clReleaseMemObject(sub_c.deviceBuffer);

    cout << "==Memory object test\n";
    auto *agent = new HeteroComputeAgent();
    const size_t dl = 1024 * 1024;
    const size_t size = dl * sizeof(float);

    cout << "==Apply buffer\n";
    auto za = agent->applyZeroCopyBuffer(size);
    printUserZBuffer(&za);
    auto zb = agent->applyZeroCopyBuffer(size);
    printUserZBuffer(&zb);
    auto zc = agent->applyZeroCopyBuffer(size);
    printUserZBuffer(&zc);
    auto zd = agent->applyZeroCopyBuffer(size);
    printUserZBuffer(&zd);

    cout << "==Gen random number\n";
    randArrayGenerator<float>(0.0, 1.0, (float *) za.hostPtr, dl);
    randArrayGenerator<float>(0.0, 1.0, (float *) zb.hostPtr, dl);

    cout << "==Release buffer\n";
    agent->releaseZeroCopyBuffer(&za);
    agent->releaseZeroCopyBuffer(&zb);
    agent->releaseZeroCopyBuffer(&zc);
    agent->releaseZeroCopyBuffer(&zd);

    return 0;
}