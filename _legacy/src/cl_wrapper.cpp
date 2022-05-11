
#include <cl_wrapper.h>

void cl_wrapper::initWrapper(cl_platform_id platform_id, cl_device_id device_id) {
    this->platform = platform_id;
    this->device = device_id;
    this->context = clCreateContext(nullptr, 1, &device_id, nullptr, nullptr, &err);
    checkError(err);
    this->memCmdQueue = createProfilingCmdQueue();
    checkError(err);
}

cl_wrapper::cl_wrapper() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    initWrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
}

cl_wrapper::cl_wrapper(int platform_num, int device_num) {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    initWrapper(clInfo.platforms[platform_num].platformId,
                clInfo.platforms[platform_num].devices[device_num].deviceId);
}

cl_wrapper::cl_wrapper(cl_platform_id platform, cl_device_id device) {
    initWrapper(platform, device);
}

cl_command_queue cl_wrapper::createCmdQueue(cl_command_queue_properties properties) {
    cl_command_queue queue = clCreateCommandQueue(context, device, properties, &err);
    checkError(err);
    this->cmdQueues.push_back(queue);
    return queue;
}

cl_command_queue cl_wrapper::createProfilingCmdQueue() {
    return createCmdQueue(CL_QUEUE_PROFILING_ENABLE);
}

cl_program cl_wrapper::createProgramWithOptions(const char **program_source, cl_uint source_len, const char *options) {
    cl_program program = clCreateProgramWithSource(context, source_len, program_source, nullptr, &err);
    checkError(err);

    //err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    err = clBuildProgram(program, 1, &this->device, options, nullptr, nullptr);
    if (err != CL_SUCCESS) {
        std::cerr << "Error " << err << " with clBuildProgram.\n";
        static const size_t LOG_SIZE = 10240;
        char log[LOG_SIZE];
        log[0] = 0;
        cl_int error = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, LOG_SIZE, log, nullptr);
        if (error == CL_INVALID_VALUE)
            std::cerr << "There was a build error, but there is insufficient space allocated to show the build logs.\n";
        else
            std::cerr << "Build error:\n" << log << "\n";
        checkError(err);
    }
    programs.push_back(program);
    return program;
}

cl_program cl_wrapper::createProgram(const char **program_source, cl_uint source_len) {
    cl_program program = clCreateProgramWithSource(context, source_len, program_source, nullptr, &err);
    checkError(err);

    //err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    err = clBuildProgram(program, 1, &this->device, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS) {
        std::cerr << "Error " << err << " with clBuildProgram.\n";
        static const size_t LOG_SIZE = 2048;
        char log[LOG_SIZE];
        log[0] = 0;
        cl_int error = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, LOG_SIZE, log, nullptr);
        if (error == CL_INVALID_VALUE)
            std::cerr << "There was a build error, but there is insufficient space allocated to show the build logs.\n";
        else
            std::cerr << "Build error:\n" << log << "\n";
        checkError(err);
    }
    programs.push_back(program);
    return program;
}

cl_kernel cl_wrapper::createKernel(const std::string &kernel_name, cl_program program) {
    cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error " << err << " create kernel " << kernel_name << std::endl;
        checkError(err);
    }
    kernels.push_back(kernel);
    return kernel;
}


cl_kernel cl_wrapper::createKernel(cl_program program, const std::string &kernel_name) {
    cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error " << err << " create kernel " << kernel_name << std::endl;
        checkError(err);
    }
    kernels.push_back(kernel);
    return kernel;
}

cl_wrapper::~cl_wrapper() {
    //Clear OpenCL
    for (auto kernel : kernels)
        clReleaseKernel(kernel);
    for (auto program : programs)
        clReleaseProgram(program);
    for (auto queue : cmdQueues)
        clReleaseCommandQueue(queue);
    //clReleaseCommandQueue(memCmdQueue);
    clReleaseContext(context);
    clReleaseDevice(this->device);
}

//static void queryCLInfo(clInfo *clinfo);
void cl_wrapper::queryCLInfo(CLInfo *clinfo) {
    cl_uint num_platforms;
    checkError(clGetPlatformIDs(0, nullptr, &num_platforms));
    clinfo->num_platforms = num_platforms;
    clinfo->platforms =
            (struct CLPlatform *) malloc(num_platforms * sizeof(struct CLPlatform));
    cl_platform_id platform_ids[num_platforms];
    checkError(clGetPlatformIDs(num_platforms, platform_ids, nullptr));
    for (int i = 0; i < num_platforms; i++) //for each platform
    {
        CLPlatform *cl_platform = &clinfo->platforms[i];
        cl_platform->platformId = platform_ids[i];
        cl_uint num_devices;
        checkError(clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_ALL, 0, nullptr, &num_devices));
        cl_platform->devices =
                (struct CLDevice *) malloc(num_devices * sizeof(struct CLDevice));
        cl_platform->num_devices = num_devices;
        cl_device_id device_ids[num_devices];
        checkError(clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_ALL, sizeof(device_ids), device_ids, NULL));
        for (int j = 0; j < num_devices; j++) //for each device
        {
            CLDevice *cl_device = &cl_platform->devices[j];
            cl_device->deviceId = device_ids[j];
            cl_ulong value = 0;
            checkError(clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &value, 0));
            cl_device->num_CUs = static_cast<int>(value);
            checkError(clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &value, 0));
            cl_device->localMemorySizeB = static_cast<int>(value);
            cl_uint dimSize;
            checkError(
                    clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &dimSize,
                                    0));
            cl_device->WGDimSize = dimSize;
            size_t work_item_sizes[dimSize];
            checkError(clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                                       sizeof(work_item_sizes), work_item_sizes, 0));
            cl_device->maxWGSizes = (size_t *) malloc(sizeof(size_t) * dimSize);
            for (size_t work_item_dim = 0; work_item_dim < dimSize; work_item_dim++)
                cl_device->maxWGSizes[work_item_dim] = work_item_sizes[work_item_dim];
            checkError(clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &value, 0));
            cl_device->maxAllocSizeB = static_cast<int>(value);
        }
    }
}

//static void printCLInfo(clInfo *clinfo); 
std::string cl_wrapper::errorMessage(cl_int error) {
    return toString(error);
}

void cl_wrapper::checkError(cl_int error) {
    if (error != CL_SUCCESS) {
        std::string message = getErrorString(error);
        std::cout << "OpenCL execution error, code " << error << " " << message << std::endl;
        throw std::runtime_error(std::string("OpenCL error, code: ") + message);
    }
}

std::string cl_wrapper::getErrorString(cl_int error) {
    switch (error) {
        // run-time and JIT compiler errors
        case 0:
            return "CL_SUCCESS";
        case -1:
            return "CL_DEVICE_NOT_FOUND";
        case -2:
            return "CL_DEVICE_NOT_AVAILABLE";
        case -3:
            return "CL_COMPILER_NOT_AVAILABLE";
        case -4:
            return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case -5:
            return "CL_OUT_OF_RESOURCES";
        case -6:
            return "CL_OUT_OF_HOST_MEMORY";
        case -7:
            return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case -8:
            return "CL_MEM_COPY_OVERLAP";
        case -9:
            return "CL_IMAGE_FORMAT_MISMATCH";
        case -10:
            return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case -11:
            return "CL_BUILD_PROGRAM_FAILURE";
        case -12:
            return "CL_MAP_FAILURE";
        case -13:
            return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case -14:
            return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case -15:
            return "CL_COMPILE_PROGRAM_FAILURE";
        case -16:
            return "CL_LINKER_NOT_AVAILABLE";
        case -17:
            return "CL_LINK_PROGRAM_FAILURE";
        case -18:
            return "CL_DEVICE_PARTITION_FAILED";
        case -19:
            return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

            // compile-time errors
        case -30:
            return "CL_INVALID_VALUE";
        case -31:
            return "CL_INVALID_DEVICE_TYPE";
        case -32:
            return "CL_INVALID_PLATFORM";
        case -33:
            return "CL_INVALID_DEVICE";
        case -34:
            return "CL_INVALID_CONTEXT";
        case -35:
            return "CL_INVALID_QUEUE_PROPERTIES";
        case -36:
            return "CL_INVALID_COMMAND_QUEUE";
        case -37:
            return "CL_INVALID_HOST_PTR";
        case -38:
            return "CL_INVALID_MEM_OBJECT";
        case -39:
            return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case -40:
            return "CL_INVALID_IMAGE_SIZE";
        case -41:
            return "CL_INVALID_SAMPLER";
        case -42:
            return "CL_INVALID_BINARY";
        case -43:
            return "CL_INVALID_BUILD_OPTIONS";
        case -44:
            return "CL_INVALID_PROGRAM";
        case -45:
            return "CL_INVALID_PROGRAM_EXECUTABLE";
        case -46:
            return "CL_INVALID_KERNEL_NAME";
        case -47:
            return "CL_INVALID_KERNEL_DEFINITION";
        case -48:
            return "CL_INVALID_KERNEL";
        case -49:
            return "CL_INVALID_ARG_INDEX";
        case -50:
            return "CL_INVALID_ARG_VALUE";
        case -51:
            return "CL_INVALID_ARG_SIZE";
        case -52:
            return "CL_INVALID_KERNEL_ARGS";
        case -53:
            return "CL_INVALID_WORK_DIMENSION";
        case -54:
            return "CL_INVALID_WORK_GROUP_SIZE";
        case -55:
            return "CL_INVALID_WORK_ITEM_SIZE";
        case -56:
            return "CL_INVALID_GLOBAL_OFFSET";
        case -57:
            return "CL_INVALID_EVENT_WAIT_LIST";
        case -58:
            return "CL_INVALID_EVENT";
        case -59:
            return "CL_INVALID_OPERATION";
        case -60:
            return "CL_INVALID_GL_OBJECT";
        case -61:
            return "CL_INVALID_BUFFER_SIZE";
        case -62:
            return "CL_INVALID_MIP_LEVEL";
        case -63:
            return "CL_INVALID_GLOBAL_WORK_SIZE";
        case -64:
            return "CL_INVALID_PROPERTY";
        case -65:
            return "CL_INVALID_IMAGE_DESCRIPTOR";
        case -66:
            return "CL_INVALID_COMPILER_OPTIONS";
        case -67:
            return "CL_INVALID_LINKER_OPTIONS";
        case -68:
            return "CL_INVALID_DEVICE_PARTITION_COUNT";

            // extension errors
        case -1000:
            return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
        case -1001:
            return "CL_PLATFORM_NOT_FOUND_KHR";
        case -1002:
            return "CL_INVALID_D3D10_DEVICE_KHR";
        case -1003:
            return "CL_INVALID_D3D10_RESOURCE_KHR";
        case -1004:
            return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
        case -1005:
            return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
        default:
            return "Unknown OpenCL error";
    }
}

double cl_wrapper::getExecutionTime(cl_event *e) {
    cl_ulong time_start, time_end;
    //clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, nullptr);
    //clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, nullptr);
    clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_QUEUED, sizeof(time_start), &time_start, nullptr);
    clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_COMPLETE, sizeof(time_end), &time_end, nullptr);
    double clExeTime = static_cast<double>(time_end - time_start) / 1000000000.0;
    //std::cout << "CL profile time: " << clExeTime << " s" << std::endl;
    //clReleaseEvent(*e);
    return clExeTime;
}




