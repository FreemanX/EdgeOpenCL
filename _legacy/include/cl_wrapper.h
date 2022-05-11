#ifndef CL_WRAPPWER_H
#define CL_WRAPPWER_H

#include <CL/cl.h>

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <map>
#include <memory>


struct CLDevice {
    cl_device_id deviceId;
    int num_CUs;
    int localMemorySizeB;
    int WGDimSize;
    size_t *maxWGSizes;
    int maxAllocSizeB;
    // more device info adds here...
};

struct CLPlatform {
    cl_platform_id platformId;
    cl_uint num_devices;
    CLDevice *devices;
};

struct CLInfo {
    cl_uint num_platforms;
    CLPlatform *platforms;
};

class cl_wrapper {
public:
    cl_context context{};
    cl_device_id device{};
    cl_platform_id platform{};
    //command queue for memory operations, not kernel execution
    cl_command_queue memCmdQueue{};

    cl_wrapper();

    cl_wrapper(cl_platform_id platform_id, cl_device_id device_id);

    cl_wrapper(int platform_num, int device_num);

    ~cl_wrapper();

    template<typename T>
    static std::string toString(T val) {
        std::ostringstream myostringstream;
        myostringstream << val;
        return myostringstream.str();
    }

    static std::string errorMessage(cl_int error);

    static void checkError(cl_int error);

    static void queryCLInfo(CLInfo *clinfo);

    static std::string getErrorString(cl_int error);

    static double getExecutionTime(cl_event *e);

    cl_command_queue createCmdQueue(cl_command_queue_properties properties);

    cl_command_queue createProfilingCmdQueue();

    cl_kernel createKernel(const std::string &kernel_name, cl_program program);

    cl_kernel createKernel(cl_program program, const std::string &kernel_name);

    cl_program createProgram(const char **program_source, cl_uint source_len);

    cl_program createProgramWithOptions(const char **program_source, cl_uint source_len, const char *options);

private:
    cl_int err{};
    std::vector<cl_command_queue> cmdQueues;
    std::vector<cl_program> programs;
    std::vector<cl_kernel> kernels;

    void initWrapper(cl_platform_id platform_id, cl_device_id device_id);
};

#endif