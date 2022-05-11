#ifndef MEMORY_MANAGER_HPP
#define MEMORY_MANAGER_HPP

#include "cl_wrapper.h"
#include <vector>
#include <map>
#include <string>


//non-zero copy method
template<typename T>
cl_mem cl_make_array(cl_wrapper *cl, T *x, size_t n);

template<typename T>
void cl_push_array(cl_wrapper *cl, cl_mem *buff, T *x, size_t n);

template<typename T>
void cl_pull_array(cl_wrapper *cl, cl_mem *buff, T *x, size_t n);

template<typename T>
cl_mem cl_make_array(cl_wrapper *cl, T *x, size_t n) {
    cl_int error;
    cl_mem buff;
    size_t size = sizeof(T) * n;
    buff = clCreateBuffer(cl->context, CL_MEM_READ_WRITE, size, nullptr, &error);
    cl_wrapper::checkError(error);

    if (x) {
        cl_push_array<T>(cl, &buff, x, n);
    }
    return buff;
}

//push array to existing cl_mem
template<typename T>
void cl_push_array(cl_wrapper *cl, cl_mem *buff, T *x, size_t n) {
    size_t size = sizeof(T) * n;
    cl_event e;
    cl_wrapper::checkError(clEnqueueWriteBuffer(cl->memCmdQueue, *buff, CL_TRUE, 0, size, x, 0, NULL, &e));
    clWaitForEvents(1, &e);
    clReleaseEvent(e);
}

template<typename T>
void cl_pull_array(cl_wrapper *cl, cl_mem *buff, T *x, size_t n) {
    cl_int error;
    cl_event e;
    size_t size = sizeof(T) * n;
    error = clEnqueueReadBuffer(cl->memCmdQueue, *buff, CL_TRUE, 0, size, x, 0, NULL, &e);
    cl_wrapper::checkError(error);
    cl_int err = clWaitForEvents(1, &e);
    clReleaseEvent(e);
    if (err != CL_SUCCESS) {
        throw std::runtime_error("wait for event on copy2host failed with " + cl_wrapper::toString(err));
    }
}

//====================Legacy code, implementation of zero copy =======================
//OpenCL managed Zero Copy memory implementation
template<typename T>
struct ZeroCopyMem {
    cl_mem deviceBuffer;
    T *hostPtr;
    size_t length; // how many elements (T)
    size_t size; // how many bytes
};

//Step 0: call init_zero_copy_region()
template<typename T>
ZeroCopyMem<T> init_zero_copy_region(cl_wrapper *cl, size_t dataLength);

//Step 1: use zeroCopyMem.hostPtr to load data @ the host side

//Step 2: call unmap_zero_copy_region()
//TODO Not sure if unmaping memory object is necessary before GPU use cl_mem
template<typename T>
void unmap_zero_copy_region(cl_wrapper *cl, ZeroCopyMem<T> *zeroCopyMem);

template<typename T>
ZeroCopyMem<T> init_zero_copy_region(cl_wrapper *cl, size_t dataLength) {
    cl_int error;
    ZeroCopyMem<T> zeroCopyMem{};
    zeroCopyMem.length = dataLength;
    zeroCopyMem.size = sizeof(T) * dataLength;
    zeroCopyMem.deviceBuffer = clCreateBuffer(cl->context,
                                              CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                              sizeof(T) * dataLength,
                                              nullptr,
                                              &error);
    cl_wrapper::checkError(error);

    //cl_event e;

    zeroCopyMem.hostPtr = (T *) clEnqueueMapBuffer(cl->memCmdQueue,
                                                   zeroCopyMem.deviceBuffer,
                                                   CL_TRUE,
                                                   CL_MAP_WRITE,
                                                   0,
                                                   sizeof(T) * dataLength,
                                                   0, NULL, NULL, &error);
    //clWaitForEvents(1, &e);
    //cl_ulong time_start, time_end;
    //clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, nullptr);
    //clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, nullptr);
    //double clExeTime = static_cast<double>(time_end - time_start)/100000000.0;
    //std::cout << "clEnqueueMapBuffer time: " << clExeTime << " s" << std::endl;

    cl_wrapper::checkError(error);
    return zeroCopyMem;
}

template<typename T>
ZeroCopyMem<T> init_zero_copy_region_byte(cl_wrapper *cl, size_t size) {
    cl_int error;
    ZeroCopyMem<T> zeroCopyMem{};
    zeroCopyMem.length = size / sizeof(T);
    zeroCopyMem.size = size;
    zeroCopyMem.deviceBuffer = clCreateBuffer(cl->context,
                                              CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
                                              size,
                                              nullptr,
                                              &error);
    cl_wrapper::checkError(error);

    map_zero_copy_region(cl, &zeroCopyMem, size);
    return zeroCopyMem;
}

template<typename T>
void map_zero_copy_region(cl_wrapper *cl, ZeroCopyMem<T> *zeroCopyMem, size_t size) {
    cl_int error;
    zeroCopyMem->hostPtr = (T *) clEnqueueMapBuffer(cl->memCmdQueue,
                                                    zeroCopyMem->deviceBuffer,
                                                    CL_TRUE,
                                                    CL_MAP_WRITE,
                                                    0,
                                                    size,
                                                    0, NULL, NULL, &error);

    cl_wrapper::checkError(error);
}

template<typename T>
void unmap_zero_copy_region(cl_wrapper *cl, ZeroCopyMem<T> *zeroCopyMem) {
    cl_int error;
    error = clEnqueueUnmapMemObject(
            cl->memCmdQueue,
            zeroCopyMem->deviceBuffer,
            (void *) zeroCopyMem->hostPtr,
            0, NULL, NULL);
    cl_wrapper::checkError(error);
}

template<typename T>
void release_zero_copy_region(ZeroCopyMem<T> *zeroCopyMem) {
    zeroCopyMem->length = 0;
    zeroCopyMem->size = 0;
    zeroCopyMem->hostPtr = nullptr;
    clReleaseMemObject(zeroCopyMem->deviceBuffer);
}

//====================Legacy code =======================

#define BUFFER_SIZE 512 * 1024 * 1024
//#define BUFFER_SIZE 8 * 1024 * 1024
#define SIZE_DIFF_THRESHOLD 8192

struct ZBuffer {
    void *hostPtr{};
    cl_mem deviceBuffer{};
    size_t totalSize{};
    size_t remainingSize{};
    size_t subBufferHead{};
    int numUsedSubBuffer = 0;
};


struct ZSubBuffer {
    void *subHostPtr;
    cl_mem deviceSubBuffer;
    ZBuffer *zBuffer;
    size_t origin;
    size_t size;
};

struct UserZBuffer {
    void *hostPtr;
    cl_mem *deviceBuffer;
    size_t size; // How many bytes
    int key;
};

void printUserZBuffer(UserZBuffer *userZBuffer);

void printZSubBuffer(ZSubBuffer *zSubBuffer);

class MemoryManager {
public:
    explicit MemoryManager(cl_wrapper *cl);

    UserZBuffer applyZeroCopyBuffer(size_t size);

    void releaseUserZBuffer(UserZBuffer *userZBuffer);

    void setDefaultBufferSize(size_t size) { DEFAULT_BUFFER_SIZE = size; }

    size_t getDefaultBufferSize() { return DEFAULT_BUFFER_SIZE; }

    ~MemoryManager();

private:
    cl_wrapper *cl;
    cl_int error{};
    cl_event clEvent{};
    cl_ulong max_mem_alloc_size{};
    cl_uint mem_base_addr_align{};
    size_t DEFAULT_BUFFER_SIZE;
    int userBufferCnt = 0;
    std::vector<ZBuffer> zBufferList;
    std::multimap<ZBuffer *, ZSubBuffer> freeZSubBufferMMap;
    std::map<int, ZSubBuffer> usedZSubBufferMap;

    void expendZBuffer();

    ZBuffer createZeroCopyBuffer(size_t size);

    ZSubBuffer createZeroCopySubBuffer(ZBuffer *zBuffer, size_t size);

    size_t convert2MemAlignedSize(size_t size);

    void insertUsedZSubBufferMap(std::pair<int, ZSubBuffer> userKeyZSubBufferPair);

    void eraseUsedZSubBufferMap(int key);

    static void assembleUserZBuffer(ZSubBuffer *zSubBuffer, UserZBuffer *userZBuffer, size_t size);

    void printMemStatus(std::string msg);

};

#endif