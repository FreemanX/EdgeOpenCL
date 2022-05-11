#pragma clang diagnostic push
#pragma ide diagnostic ignored "readability-convert-member-functions-to-static"
#pragma ide diagnostic ignored "OCDFAInspection"
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include "CLRT.h"
#include <string>
#include <utility>
#include <cstdio>
#include <dlfcn.h>
#include "utils.h"

void CLRT::
initRT(cl_platform_id Platform, cl_device_id Device) {
	this->platform_id_ = Platform;
	this->device_id_ = Device;
	this->context_ = clCreateContext(nullptr, 1, &device_id_, nullptr, nullptr, &err, oclLibHandle_);
	CLRT_ERR(err);
	this->memory_queue_ = clCreateCommandQueue(context_, device_id_, 0, &err, oclLibHandle_);
	CLRT_ERR(err);
	initSubDevices();
}

void CLRT::queryCLInfo(InternalInfo_CL_ *clinfo) const {
	cl_uint num_platforms;
	CLRT_ERR(clGetPlatformIDs(0, nullptr, &num_platforms, oclLibHandle_));
	clinfo->num_platforms = num_platforms;
	clinfo->platforms =
			(InternalPlatform_CL_ *)malloc(num_platforms * sizeof(InternalPlatform_CL_));
	cl_platform_id platform_ids[num_platforms];
	CLRT_ERR(clGetPlatformIDs(num_platforms, platform_ids, nullptr, oclLibHandle_));
	for (int i = 0; i < num_platforms; i++) { //for each platform
		InternalPlatform_CL_ *cl_platform = &clinfo->platforms[i];
		cl_platform->platformId = platform_ids[i];
		cl_uint num_devices;
		CLRT_ERR(clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_ALL, 0, nullptr, &num_devices, oclLibHandle_));
		cl_platform->devices = (InternalDevice_CL_ *)malloc(num_devices * sizeof(InternalDevice_CL_));
		cl_platform->num_devices = num_devices;
		cl_device_id device_ids[num_devices];
		CLRT_ERR(clGetDeviceIDs(platform_ids[i],
														CL_DEVICE_TYPE_ALL,
														sizeof(device_ids),
														device_ids,
														nullptr,
														oclLibHandle_));
		for (int j = 0; j < num_devices; j++) { //for each device
			InternalDevice_CL_ *cl_device = &cl_platform->devices[j];
			cl_device->deviceId = device_ids[j];
			cl_ulong value = 0;
			CLRT_ERR(clGetDeviceInfo(cl_device->deviceId,
															 CL_DEVICE_MAX_COMPUTE_UNITS,
															 sizeof(cl_ulong),
															 &value,
															 nullptr,
															 oclLibHandle_));
			cl_device->num_CUs = static_cast<int>(value);
			CLRT_ERR(clGetDeviceInfo(cl_device->deviceId,
															 CL_DEVICE_LOCAL_MEM_SIZE,
															 sizeof(cl_ulong),
															 &value,
															 nullptr,
															 oclLibHandle_));
			cl_device->localMemorySizeB = static_cast<int>(value);
			cl_uint dimSize;
			CLRT_ERR(
					clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint),
													&dimSize, nullptr, oclLibHandle_));
			cl_device->WGDimSize = dimSize;
			size_t work_item_sizes[dimSize];
			CLRT_ERR(clGetDeviceInfo(cl_device->deviceId, CL_DEVICE_MAX_WORK_ITEM_SIZES,
															 sizeof(work_item_sizes), work_item_sizes, nullptr, oclLibHandle_));
			cl_device->maxWGSizes = (size_t *)malloc(sizeof(size_t) * dimSize);
			for (size_t work_item_dim = 0; work_item_dim < dimSize; work_item_dim++)
				cl_device->maxWGSizes[work_item_dim] = work_item_sizes[work_item_dim];
			CLRT_ERR(clGetDeviceInfo(cl_device->deviceId,
															 CL_DEVICE_MAX_MEM_ALLOC_SIZE,
															 sizeof(cl_ulong),
															 &value,
															 nullptr,
															 oclLibHandle_));
			cl_device->maxAllocSizeB = static_cast<int>(value);
			cl_uint sizeBufferInt;
			clGetDeviceInfo(cl_device->deviceId,
											CL_DEVICE_PARTITION_MAX_SUB_DEVICES,
											sizeof(sizeBufferInt),
											&sizeBufferInt,
											nullptr, oclLibHandle_);
			cl_device->maxSubDevices = sizeBufferInt;
		}
	}
}

CLRT::CLRT() {
	this->oclLibPath_ = DEFAULT_GPU_OCL_PATH;
	open_opencl_so(oclLibPath_.c_str(), &oclLibHandle_);
	queryCLInfo(&clInfo);
	initRT(clInfo.platforms[0].platformId, clInfo.platforms->devices[0].deviceId);
}

CLRT::CLRT(std::string OclLibPath) {
	this->oclLibPath_ = std::move(OclLibPath);
	open_opencl_so(oclLibPath_.c_str(), &oclLibHandle_);
	queryCLInfo(&clInfo);
	initRT(clInfo.platforms[0].platformId, clInfo.platforms->devices[0].deviceId);
}

std::string CLRT::
CLRTErrMsg(cl_int error) {
	switch (error) {
		// run-time and JIT compiler errors
		case 0:return "CL_SUCCESS";
		case -1:return "CL_DEVICE_NOT_FOUND";
		case -2:return "CL_DEVICE_NOT_AVAILABLE";
		case -3:return "CL_COMPILER_NOT_AVAILABLE";
		case -4:return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case -5:return "CL_OUT_OF_RESOURCES";
		case -6:return "CL_OUT_OF_HOST_MEMORY";
		case -7:return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case -8:return "CL_MEM_COPY_OVERLAP";
		case -9:return "CL_IMAGE_FORMAT_MISMATCH";
		case -10:return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case -11:return "CL_BUILD_PROGRAM_FAILURE";
		case -12:return "CL_MAP_FAILURE";
		case -13:return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case -14:return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case -15:return "CL_COMPILE_PROGRAM_FAILURE";
		case -16:return "CL_LINKER_NOT_AVAILABLE";
		case -17:return "CL_LINK_PROGRAM_FAILURE";
		case -18:return "CL_DEVICE_PARTITION_FAILED";
		case -19:return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
			// compile-time errors
		case -30:return "CL_INVALID_VALUE";
		case -31:return "CL_INVALID_DEVICE_TYPE";
		case -32:return "CL_INVALID_PLATFORM";
		case -33:return "CL_INVALID_DEVICE";
		case -34:return "CL_INVALID_CONTEXT";
		case -35:return "CL_INVALID_QUEUE_PROPERTIES";
		case -36:return "CL_INVALID_COMMAND_QUEUE";
		case -37:return "CL_INVALID_HOST_PTR";
		case -38:return "CL_INVALID_MEM_OBJECT";
		case -39:return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case -40:return "CL_INVALID_IMAGE_SIZE";
		case -41:return "CL_INVALID_SAMPLER";
		case -42:return "CL_INVALID_BINARY";
		case -43:return "CL_INVALID_BUILD_OPTIONS";
		case -44:return "CL_INVALID_PROGRAM";
		case -45:return "CL_INVALID_PROGRAM_EXECUTABLE";
		case -46:return "CL_INVALID_KERNEL_NAME";
		case -47:return "CL_INVALID_KERNEL_DEFINITION";
		case -48:return "CL_INVALID_KERNEL";
		case -49:return "CL_INVALID_ARG_INDEX";
		case -50:return "CL_INVALID_ARG_VALUE";
		case -51:return "CL_INVALID_ARG_SIZE";
		case -52:return "CL_INVALID_KERNEL_ARGS";
		case -53:return "CL_INVALID_WORK_DIMENSION";
		case -54:return "CL_INVALID_WORK_GROUP_SIZE";
		case -55:return "CL_INVALID_WORK_ITEM_SIZE";
		case -56:return "CL_INVALID_GLOBAL_OFFSET";
		case -57:return "CL_INVALID_EVENT_WAIT_LIST";
		case -58:return "CL_INVALID_EVENT";
		case -59:return "CL_INVALID_OPERATION";
		case -60:return "CL_INVALID_GL_OBJECT";
		case -61:return "CL_INVALID_BUFFER_SIZE";
		case -62:return "CL_INVALID_MIP_LEVEL";
		case -63:return "CL_INVALID_GLOBAL_WORK_SIZE";
		case -64:return "CL_INVALID_PROPERTY";
		case -65:return "CL_INVALID_IMAGE_DESCRIPTOR";
		case -66:return "CL_INVALID_COMPILER_OPTIONS";
		case -67:return "CL_INVALID_LINKER_OPTIONS";
		case -68:return "CL_INVALID_DEVICE_PARTITION_COUNT";
			// extension errors
		case -1000:return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
		case -1001:return "CL_PLATFORM_NOT_FOUND_KHR";
		case -1002:return "CL_INVALID_D3D10_DEVICE_KHR";
		case -1003:return "CL_INVALID_D3D10_RESOURCE_KHR";
		case -1004:return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
		case -1005:return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
		default:return "Unknown OpenCL error";
	}
}

void CLRT::
CLRTChkErr(cl_int error) {
	if (error != CL_SUCCESS) {
		std::string message = CLRTErrMsg(error);
		std::cerr << "OpenCL execution error, code " << error << " " << message << std::endl;
		throw std::runtime_error(std::string("OpenCL error ") + message);
	}
}

double CLRT::
getExeTime(cl_event const *e) const {
	if (e == nullptr) flushed_printf("Can't get event time from null event!\n");
	cl_ulong time_start, time_end;
	clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, nullptr, oclLibHandle_);
	clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, nullptr, oclLibHandle_);
//  clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_QUEUED, sizeof(time_start), &time_start, nullptr, handle);
//  clGetEventProfilingInfo(*e, CL_PROFILING_COMMAND_COMPLETE, sizeof(time_end), &time_end, nullptr, handle);
	double clExeTime = static_cast<double>(time_end - time_start) / 1000000000.0;
//  std::cout << "Time start: " << time_start << ", time end: " << time_end << std::endl;
//  std::cout << "CL profile time: " << clExeTime << " s" << std::endl;
	//clReleaseEvent(*e);
	return clExeTime;
}

void CLRT::
printCLInfo() const {
	cl_platform_id platform = this->platform_id_;
	cl_device_id device = this->device_id_;
	char buffer[1024];
	clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(buffer), buffer, nullptr, oclLibHandle_);
	std::cout << "Platform name: " << buffer << std::endl;
	clGetPlatformInfo(platform, CL_PLATFORM_PROFILE, sizeof(buffer), buffer, nullptr, oclLibHandle_);
	std::cout << "CL_PLATFORM_PROFILE: " << buffer << std::endl;
	clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(buffer), buffer, nullptr, oclLibHandle_);
	std::cout << "CL_PLATFORM_VERSION: " << buffer << std::endl;
	clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(buffer), buffer, nullptr, oclLibHandle_);
	std::cout << "CL_PLATFORM_VENDOR: " << buffer << std::endl;
	clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS, sizeof(buffer), buffer, nullptr, oclLibHandle_);
	std::cout << "CL_PLATFORM_EXTENSIONS: " << buffer << std::endl;

	char device_string[1024];

	// CL_DEVICE_NAME
	clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_string), &device_string, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_NAME: \t\t\t%s\n", device_string);

	// CL_DEVICE_VENDOR
	clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(device_string), &device_string, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_VENDOR: \t\t\t%s\n", device_string);

	// CL_DRIVER_VERSION
	clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(device_string), &device_string, nullptr, oclLibHandle_);
	printf("  CL_DRIVER_VERSION: \t\t\t%s\n", device_string);

	// CL_DRIVER_VERSION
	clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_EXTENSIONS: \t\t%s\n", device_string);

	// CL_DEVICE_INFO
	cl_device_type type;
	clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(type), &type, nullptr, oclLibHandle_);
	if (type & CL_DEVICE_TYPE_CPU)
		printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_CPU");
	if (type & CL_DEVICE_TYPE_GPU)
		printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_GPU");
	if (type & CL_DEVICE_TYPE_ACCELERATOR)
		printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_ACCELERATOR");
	if (type & CL_DEVICE_TYPE_DEFAULT)
		printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_DEFAULT");

	// CL_DEVICE_MAX_COMPUTE_UNITS
	cl_uint compute_units;
	clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_MAX_COMPUTE_UNITS:\t\t%u\n", compute_units);

	// CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
	size_t workitem_dims;
	clGetDeviceInfo(device,
									CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
									sizeof(workitem_dims),
									&workitem_dims,
									nullptr,
									oclLibHandle_);
	printf("  CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%zu\n", workitem_dims);

	// CL_DEVICE_MAX_WORK_ITEM_SIZES
	size_t workitem_size[3];
	clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_MAX_WORK_ITEM_SIZES:\t%zu / %zu / %zu \n", workitem_size[0], workitem_size[1],
				 workitem_size[2]);

	// CL_DEVICE_MAX_WORK_GROUP_SIZE
	size_t workgroup_size;
	clGetDeviceInfo(device,
									CL_DEVICE_MAX_WORK_GROUP_SIZE,
									sizeof(workgroup_size),
									&workgroup_size,
									nullptr,
									oclLibHandle_);
	printf("  CL_DEVICE_MAX_WORK_GROUP_SIZE:\t%zu\n", workgroup_size);

	// CL_DEVICE_MAX_CLOCK_FREQUENCY
	cl_uint clock_frequency;
	clGetDeviceInfo(device,
									CL_DEVICE_MAX_CLOCK_FREQUENCY,
									sizeof(clock_frequency),
									&clock_frequency,
									nullptr,
									oclLibHandle_);
	printf("  CL_DEVICE_MAX_CLOCK_FREQUENCY:\t%u MHz\n", clock_frequency);

	// CL_DEVICE_ADDRESS_BITS
	cl_uint addr_bits;
	clGetDeviceInfo(device, CL_DEVICE_ADDRESS_BITS, sizeof(addr_bits), &addr_bits, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_ADDRESS_BITS:\t\t%u\n", addr_bits);

	// CL_DEVICE_MAX_MEM_ALLOC_SIZE
	cl_ulong max_mem_alloc_size;
	clGetDeviceInfo(device,
									CL_DEVICE_MAX_MEM_ALLOC_SIZE,
									sizeof(max_mem_alloc_size),
									&max_mem_alloc_size,
									nullptr,
									oclLibHandle_);
	printf("  CL_DEVICE_MAX_MEM_ALLOC_SIZE:\t\t%u MByte\n", (unsigned int)(max_mem_alloc_size / (1024 * 1024)));

	// CL_DEVICE_GLOBAL_MEM_SIZE
	cl_ulong mem_size;
	clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_GLOBAL_MEM_SIZE:\t\t%u MByte\n", (unsigned int)(mem_size / (1024 * 1024)));

	// CL_DEVICE_ERROR_CORRECTION_SUPPORT
	cl_bool error_correction_support;
	clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support),
									&error_correction_support, nullptr);
	printf("  CL_DEVICE_ERROR_CORRECTION_SUPPORT:\t%s\n", error_correction_support == CL_TRUE ? "yes" : "no");

	// CL_DEVICE_LOCAL_MEM_TYPE
	cl_device_local_mem_type local_mem_type;
	clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_LOCAL_MEM_TYPE:\t\t%s\n", local_mem_type == 1 ? "local" : "global");

	// CL_DEVICE_LOCAL_MEM_SIZE
	clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_LOCAL_MEM_SIZE:\t\t%u KByte\n", (unsigned int)(mem_size / 1024));

	// CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
	clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE:\t%u KByte\n", (unsigned int)(mem_size / 1024));

	// CL_DEVICE_QUEUE_PROPERTIES
	cl_command_queue_properties queue_properties;
	clGetDeviceInfo(device,
									CL_DEVICE_QUEUE_PROPERTIES,
									sizeof(queue_properties),
									&queue_properties,
									nullptr,
									oclLibHandle_);
	if (queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)
		printf("  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE");
	if (queue_properties & CL_QUEUE_PROFILING_ENABLE)
		printf("  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_PROFILING_ENABLE");

	// CL_DEVICE_IMAGE_SUPPORT
	cl_bool image_support;
	clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(image_support), &image_support, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_IMAGE_SUPPORT:\t\t%u\n", image_support);

	// CL_DEVICE_MAX_READ_IMAGE_ARGS
	cl_uint max_read_image_args;
	clGetDeviceInfo(device,
									CL_DEVICE_MAX_READ_IMAGE_ARGS,
									sizeof(max_read_image_args),
									&max_read_image_args,
									nullptr,
									oclLibHandle_);
	printf("  CL_DEVICE_MAX_READ_IMAGE_ARGS:\t%u\n", max_read_image_args);

	// CL_DEVICE_MAX_WRITE_IMAGE_ARGS
	cl_uint max_write_image_args;
	clGetDeviceInfo(device, CL_DEVICE_MAX_WRITE_IMAGE_ARGS, sizeof(max_write_image_args), &max_write_image_args,
									nullptr);
	printf("  CL_DEVICE_MAX_WRITE_IMAGE_ARGS:\t%u\n", max_write_image_args);

	// CL_DEVICE_IMAGE2D_MAX_WIDTH, CL_DEVICE_IMAGE2D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_WIDTH, CL_DEVICE_IMAGE3D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_DEPTH
	size_t szMaxDims[5];
	printf("\n  CL_DEVICE_IMAGE <dim>");
	clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(size_t), &szMaxDims[0], nullptr, oclLibHandle_);
	printf("\t\t\t2D_MAX_WIDTH\t %zu\n", szMaxDims[0]);
	clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[1], nullptr, oclLibHandle_);
	printf("\t\t\t\t\t2D_MAX_HEIGHT\t %zu\n", szMaxDims[1]);
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(size_t), &szMaxDims[2], nullptr, oclLibHandle_);
	printf("\t\t\t\t\t3D_MAX_WIDTH\t %zu\n", szMaxDims[2]);
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[3], nullptr, oclLibHandle_);
	printf("\t\t\t\t\t3D_MAX_HEIGHT\t %zu\n", szMaxDims[3]);
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(size_t), &szMaxDims[4], nullptr, oclLibHandle_);
	printf("\t\t\t\t\t3D_MAX_DEPTH\t %zu\n", szMaxDims[4]);

	// CL_DEVICE_PREFERRED_VECTOR_WIDTH_<type>
	printf("  CL_DEVICE_PREFERRED_VECTOR_WIDTH_<t>\t");
	cl_uint vec_width[6];
	clGetDeviceInfo(device,
									CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,
									sizeof(cl_uint),
									&vec_width[0],
									nullptr,
									oclLibHandle_);
	clGetDeviceInfo(device,
									CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,
									sizeof(cl_uint),
									&vec_width[1],
									nullptr,
									oclLibHandle_);
	clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(cl_uint), &vec_width[2], nullptr, oclLibHandle_);
	clGetDeviceInfo(device,
									CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,
									sizeof(cl_uint),
									&vec_width[3],
									nullptr,
									oclLibHandle_);
	clGetDeviceInfo(device,
									CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,
									sizeof(cl_uint),
									&vec_width[4],
									nullptr,
									oclLibHandle_);
	clGetDeviceInfo(device,
									CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,
									sizeof(cl_uint),
									&vec_width[5],
									nullptr,
									oclLibHandle_);
	printf("CHAR %u, SHORT %u, INT %u, FLOAT %u, DOUBLE %u\n",
				 vec_width[0], vec_width[1], vec_width[2], vec_width[3], vec_width[4]);

	printf("\n");
	cl_uint sizeBufferInt;
	clGetDeviceInfo(device,
									CL_DEVICE_PARTITION_MAX_SUB_DEVICES,
									sizeof(sizeBufferInt),
									&sizeBufferInt,
									nullptr,
									oclLibHandle_);
	std::cout << "  CL_DEVICE_PARTITION_MAX_SUB_DEVICES:\t" << sizeBufferInt << std::endl;

	cl_device_partition_property properties[5];
	for (long p : properties) { p = 0; }
	clGetDeviceInfo(device, CL_DEVICE_PARTITION_PROPERTIES, sizeof(properties), properties, nullptr, oclLibHandle_);
	std::cout << "  CL_DEVICE_PARTITION_PROPERTIES: \n";
	for (long p : properties) {
		switch (p) {
			case CL_DEVICE_PARTITION_EQUALLY: std::cout << "\t\t\t\t\tCL_DEVICE_PARTITION_EQUALLY " << p << "\n";
			case CL_DEVICE_PARTITION_BY_COUNTS: std::cout << "\t\t\t\t\tCL_DEVICE_PARTITION_BY_COUNTS " << p << "\n";
			case CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN:
				std::cout << "\t\t\t\t\tCL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN " << p << "\n";
			default: std::cout << "\t\t\t\t\tUnknown: " << p << " " << "\n";
		}
	}
	printf("\n");

	cl_uint addrAlign;
	clGetDeviceInfo(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(addrAlign), &addrAlign, nullptr, oclLibHandle_);
	std::cout << "  CL_DEVICE_MEM_BASE_ADDR_ALIGN:\t" << addrAlign << std::endl;

	cl_device_svm_capabilities svm;
	clGetDeviceInfo(device, CL_DEVICE_SVM_CAPABILITIES, sizeof(svm), &svm, nullptr, oclLibHandle_);
	printf("  CL_DEVICE_SVM_CAPABILITIES: ");
	printf("\t\tCL_DEVICE_SVM_FINE_GRAIN_SYSTEM: %lu\n",
				 (svm & CL_DEVICE_SVM_FINE_GRAIN_SYSTEM));
	printf("\t\t\t\t\tCL_DEVICE_SVM_FINE_GRAIN_BUFFER %lu\n",
				 (svm & CL_DEVICE_SVM_FINE_GRAIN_BUFFER));
	printf("\t\t\t\t\tCL_DEVICE_SVM_COARSE_GRAIN_BUFFER %lu\n",
				 (svm & CL_DEVICE_SVM_COARSE_GRAIN_BUFFER));
	printf("\n\n\n");

}

cl_program CLRT::
createProgramFromSource_cl(cl_context Context,
													 cl_device_id *Devices,
													 unsigned int NumDevices,
													 const char **Source,
													 const char *Options) {
	const size_t source_len = strlen(*Source);
	cl_program program = clCreateProgramWithSource(Context, 1, Source, &source_len, &err, oclLibHandle_);
	CLRT_ERR(err);
	err = clBuildProgram(program, NumDevices, Devices, Options, nullptr, nullptr, oclLibHandle_);
	for (int di = 0; di < NumDevices; ++di) {
		if (err != CL_SUCCESS) {
			std::cerr << "Error " << err << " with clBuildProgram.\n";
			static const size_t LOG_SIZE = 1048576;
//	  char log[LOG_SIZE];
			char *log = (char *)malloc(LOG_SIZE);
			log[0] = 0;
			cl_int error = clGetProgramBuildInfo(program, Devices[di], CL_PROGRAM_BUILD_LOG,
																					 LOG_SIZE, log, nullptr, oclLibHandle_);
			if (error == CL_INVALID_VALUE)
				std::cerr << "There was a build error, but there is insufficient space allocated to show the build logs.\n";
			else
				flushed_printf("Build error:\n%s\nKernel Source:\n%s\n", log, *Source);
			free(log);
			CLRT_ERR(err);
		}
	}
	return program;
}

CLRT_Program CLRT::
createProgramFromSource(const char **ProgramSource) {
	auto clrt_program = new InternalProgram_;
	clrt_program->program_ = createProgramFromSource_cl(context_, &device_id_, 1, ProgramSource, nullptr);
	programs[clrt_program->getTrackNum()] = clrt_program;
	return clrt_program;
}

CLRT_Program CLRT::
createProgramFromSourceWithOptions(const char **ProgramSource, const char *Options) {
	auto clrt_program = new InternalProgram_;
	clrt_program->program_ = createProgramFromSource_cl(context_,
																											&device_id_,
																											1,
																											ProgramSource,
																											Options);
	programs[clrt_program->getTrackNum()] = clrt_program;
	return clrt_program;
}

cl_kernel CLRT::
createKernel_cl(cl_program Program, std::string &Name) {
	cl_kernel kernel = clCreateKernel(Program, Name.c_str(), &err, oclLibHandle_);
	if (err != CL_SUCCESS)
		printf("Create Kernel %s failed with error code %d\n", Name.c_str(), err);
	CLRT_ERR(err);
	return kernel;
}

CLRT_Kernel CLRT::
createKernel(CLRT_Program &Program, std::string Name, unsigned int NumArgs) {
	auto clrt_kernel = new InternalKernel_;
	clrt_kernel->kernelName = std::move(Name);
	clrt_kernel->numArgs = NumArgs;
	clrt_kernel->kernel_ = createKernel_cl(Program->program_, clrt_kernel->kernelName);
	clrt_kernel->oclLibHandle = this->oclLibHandle_;
	clrt_kernel->argSizeList = static_cast<size_t *>(calloc(NumArgs, sizeof(size_t)));
	clrt_kernel->argPtrList = static_cast<void **>(calloc(NumArgs, sizeof(void *)));
	clrt_kernel->argTypeList = static_cast<BufferType *>(calloc(NumArgs, sizeof(BufferType)));
	kernels[clrt_kernel->getTrackNum()] = clrt_kernel;
	return clrt_kernel;
}

void CLRT::
setKernelArgWithLists(cl_kernel CLKernel,
											const BufferType *TypeList,
											void **PtrList,
											size_t *SizeList,
											uint NumArgs) const {
	for (int i = 0; i < NumArgs; ++i) {
		if (!PtrList[i]) {
			flushed_printf("Pointer is null for %dth arg!\n", i);
			exit(-1);
		}
		// print debug info
//#define PRINT_ARGS
		switch (TypeList[i]) {
			case BUF: {
#ifdef PRINT_ARGS
				flushed_printf("Arg %d, BUF, size: %d\n", i, SizeList[i]);
#endif
				auto buffer = (CLRT_Buffer)PtrList[i];
				CLRT_ERR(
						clSetKernelArg(CLKernel,
													 i,
													 sizeof(cl_mem),
													 buffer->getCL_memPtr(), oclLibHandle_));
				break;
			}
/*
	  case IMG:
#ifdef PRINT_ARGS
		flushed_printf("Arg %d, IMG, size: %d\n", i, SizeList[i]);
#endif
		CLRT_ERR(clSetKernelArg(CLRTKernel->kernel_, i, sizeof(cl_mem), PtrList[i]));
		break;
*/
			case Z_BUF: {
#ifdef PRINT_ARGS
				flushed_printf("Arg %d, Z_BUF, size: %d\n", i, SizeList[i]);
#endif
				auto zBuffer = (CLRT_zBuffer)PtrList[i];
				CLRT_ERR(
						clSetKernelArg(CLKernel,
													 i,
													 sizeof(cl_mem),
													 zBuffer->getCL_memPtr(), oclLibHandle_));
				break;
			}
			case PRIMITIVE: {
#ifdef PRINT_ARGS
				flushed_printf("Arg %d, PRIMITIVE, size: %d\n", i, SizeList[i]);
#endif
				CLRT_ERR(clSetKernelArg(CLKernel, i, SizeList[i], PtrList[i], oclLibHandle_));
				break;
			}
//	  case RELEASED: {
//		flushed_printf("Arg %d, Error!, released Buffer, size: %d\n", i, SizeList[i]);
//		exit(-1);
//	  }
			default: {
				flushed_printf("ARG %d, Error! Unknown Type! size: %d\n", i, SizeList[i]);
				exit(-1);
			}
		}
	}
}

void CLRT::
setKernelArgWithLists(CLRT_Kernel CLRTKernel,
											const BufferType *TypeList,
											void **PtrList,
											size_t *SizeList,
											uint NumArgs) const {
	setKernelArgWithLists(CLRTKernel->kernel_, TypeList, PtrList, SizeList, NumArgs);
}

void CLRT::
setKernelArgAtIndex(CLRT_Kernel Kernel,
										uint Index,
										BufferType ArgType,
										size_t ArgSize,
										void *ArgPtr) const {
	if (!ArgPtr) {
		std::cerr << "Setting null argument to " << Index << "th arg of " << Kernel->kernelName << "\n";
		exit(-1);
	}
	if (ArgType == BUF) {
		auto buffer = (CLRT_Buffer)ArgPtr;
		if (buffer->SIGNATURE == BUF) {
			CLRT_ERR(clSetKernelArg(Kernel->kernel_, Index, sizeof(cl_mem), &(buffer->mem), oclLibHandle_));
			Kernel->argTypeList[Index] = BUF;
		}
	} else if (ArgType == Z_BUF) {
		auto zbuffer = (CLRT_zBuffer)ArgPtr;
		if (zbuffer->SIGNATURE == Z_BUF) {
			CLRT_ERR(clSetKernelArg(Kernel->kernel_, Index, sizeof(cl_mem), &(zbuffer->mem), oclLibHandle_));
			Kernel->argTypeList[Index] = Z_BUF;
		}
	} else {
		CLRT_ERR(clSetKernelArg(Kernel->kernel_, Index, ArgSize, ArgPtr, oclLibHandle_));
		Kernel->argTypeList[Index] = PRIMITIVE;
	}
	Kernel->argSizeList[Index] = ArgSize;
	Kernel->argPtrList[Index] = ArgPtr;
}

void CLRT::
setKernelArgUsedByMacro_CLRT(CLRT_Kernel Kernel, size_t ArgSize, void *ArgPtr) {
	if (ArgPtr) { // ArgPtr can be NULL to skip set certain args with macro
		uint index = Kernel->argIndex % Kernel->numArgs;
		void *oclLibHandle = Kernel->oclLibHandle;
		if (ArgSize == sizeof(InternalBuffer_)) {
			auto buffer = (CLRT_Buffer)ArgPtr;
			if (buffer->SIGNATURE == BUF) {
				CLRT_ERR(clSetKernelArg(Kernel->kernel_,
																index,
																sizeof(cl_mem),
																&(buffer->mem), oclLibHandle));
				Kernel->argSizeList[index] = buffer->sizeByte;
//		Kernel->argPtrList[index] = &(buffer->mem);
				Kernel->argPtrList[index] = ArgPtr;
				Kernel->argTypeList[index] = BUF;
			}
		} else if (ArgSize == sizeof(InternalZBuffer_)) {
			auto zbuffer = (CLRT_zBuffer)ArgPtr;
			if (zbuffer->SIGNATURE == Z_BUF) {
				CLRT_ERR(clSetKernelArg(Kernel->kernel_,
																index,
																sizeof(cl_mem),
																&(zbuffer->mem), Kernel->oclLibHandle));
				Kernel->argSizeList[index] = zbuffer->sizeByte;
				Kernel->argPtrList[index] = ArgPtr;
//		Kernel->argPtrList[index] = &(zbuffer->mem);
				Kernel->argTypeList[index] = Z_BUF;
			}
		} else {
			CLRT_ERR(clSetKernelArg(Kernel->kernel_,
															index,
															ArgSize,
															ArgPtr, Kernel->oclLibHandle));
			Kernel->argSizeList[index] = ArgSize;
			Kernel->argPtrList[index] = ArgPtr;
			Kernel->argTypeList[index] = PRIMITIVE;
		}
	}
	Kernel->argIndex++;
}

void CLRT::
execKernel(cl_command_queue ExeQ,
					 CLRT_Kernel Kernel,
					 unsigned int WorkDim,
					 const size_t *Offset,
					 const size_t *GlobalSize,
					 const size_t *LocalSize,
					 unsigned int NumEventsWait,
					 const cl_event *EventWaitList,
					 cl_event *event) const {
	CLRT_ERR(clEnqueueNDRangeKernel(ExeQ,
																	Kernel->kernel_,
																	WorkDim,
																	Offset,
																	GlobalSize,
																	LocalSize,
																	NumEventsWait,
																	EventWaitList,
																	event,
																	oclLibHandle_));
}

CLRT_Buffer CLRT::
createBuffer(size_t nByte, void *HostPtr) {
	auto buffer = new InternalBuffer_(nByte);
	buffer->SIGNATURE = BUF;
	buffer->sizeByte = nByte;
	buffer->mem = clCreateBuffer(context_, CL_MEM_READ_WRITE, nByte, nullptr, &err, oclLibHandle_);
	CLRT_ERR(err);
	if (HostPtr) pushBuffer(buffer, 0, nByte, HostPtr);
	buffers[buffer->getTrackNum()] = buffer;
	return buffer;
}

CLRT_zBuffer CLRT::
createZBuffer(size_t nByte, void *HostPtr) {
	auto zbuffer = new InternalZBuffer_(nByte);
	zbuffer->SIGNATURE = Z_BUF;
	zbuffer->sizeByte = nByte;
	zbuffer->mem =
			clCreateBuffer(context_, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, nByte, nullptr, &err, oclLibHandle_);
	CLRT_ERR(err);
	mapZBuffer(zbuffer);
	if (HostPtr) memcpy(zbuffer->hostPtr, HostPtr, nByte);
	zbuffers[zbuffer->getTrackNum()] = zbuffer;
	return zbuffer;
}

void *CLRT::
mapZBuffer(CLRT_zBuffer ZBuffer) {
	if (ZBuffer->hostPtr) return ZBuffer->hostPtr;
	else {
		ZBuffer->hostPtr = clEnqueueMapBuffer(memory_queue_, ZBuffer->mem, CL_TRUE,
																					CL_MAP_WRITE | CL_MAP_READ, 0, ZBuffer->sizeByte,
																					0, nullptr, nullptr, &err, oclLibHandle_);
		clFinish(memory_queue_, oclLibHandle_);
		CLRT_ERR(err);
		return ZBuffer->hostPtr;
	}
}

void CLRT::
unmapZBuffer(CLRT_zBuffer ZBuffer) {
	if (ZBuffer->hostPtr) {
		CLRT_ERR(
				clEnqueueUnmapMemObject(memory_queue_, ZBuffer->mem, ZBuffer->hostPtr, 0, nullptr, nullptr, oclLibHandle_)
		);
		clFinish(memory_queue_, oclLibHandle_);
		ZBuffer->hostPtr = nullptr;
	}
}

void CLRT::
pushBuffer(CLRT_Buffer Buffer, size_t Offset, size_t nByte, void *HostPtr) {
	CLRT_ERR(clEnqueueWriteBuffer(memory_queue_, Buffer->mem, CL_TRUE,
																Offset, nByte, HostPtr,
																0, nullptr, nullptr, oclLibHandle_));
	clFinish(memory_queue_, oclLibHandle_);
}

void CLRT::
pullBuffer(CLRT_Buffer Buffer, size_t Offset, size_t nByte, void *HostPtr) {
	CLRT_ERR(clEnqueueReadBuffer(memory_queue_, Buffer->mem, CL_TRUE,
															 Offset, nByte, HostPtr,
															 0, nullptr, nullptr, oclLibHandle_));
	clFinish(memory_queue_, oclLibHandle_);
}

void CLRT::
releaseCLRT_Program(CLRT_Program Program) {
	if (programs.count(Program->getTrackNum()) > 0) {
		clReleaseProgram(Program->program_, oclLibHandle_);
		programs.erase(Program->getTrackNum());
		delete Program;
	}
}

void CLRT::
releaseCLRT_Kernel(CLRT_Kernel Kernel) {
	if (kernels.count(Kernel->getTrackNum()) > 0) {
		Kernel->numArgs = 0;
		delete Kernel->argSizeList;
		delete Kernel->argPtrList;
		delete Kernel->argTypeList;
		Kernel->argSizeList = nullptr;
		Kernel->argPtrList = nullptr;
		Kernel->argTypeList = nullptr;
		kernels.erase(Kernel->getTrackNum());
		clReleaseKernel(Kernel->kernel_, oclLibHandle_);
		delete Kernel;
	}
}

void CLRT::
releaseCLRT_Buffer(CLRT_Buffer Buffer) {
	if (buffers.count(Buffer->getTrackNum()) > 0) {
		clReleaseMemObject(Buffer->mem, oclLibHandle_);
		buffers.erase(Buffer->getTrackNum());
		delete Buffer;
	}
}

void CLRT::
releaseCLRT_zBuffer(CLRT_zBuffer zBuffer) {
	if (zbuffers.count(zBuffer->getTrackNum()) > 0) {
		clReleaseMemObject(zBuffer->mem, oclLibHandle_);
		zbuffers.erase(zBuffer->getTrackNum());
//	debug_printf(__func__, "%d\n", zBuffer->getTrackNum());
		delete zBuffer;
	}
}

CLRT::~CLRT() {
	for (std::vector<cl_device_id> *subDevices : sub_deviceCollection) {
		for (cl_device_id d : *subDevices) clReleaseDevice(d, oclLibHandle_);
		delete subDevices;
	}

	std::vector<int> keys;
//  flushed_printf("buffers...\n");
	for (auto &p : buffers) keys.push_back(p.first);
	for (auto &k : keys) releaseCLRT_Buffer(buffers[k]);
	keys.clear();


//  flushed_printf("zbuffers...\n");
	for (auto &p : zbuffers) keys.push_back(p.first);
	for (auto &k : keys) releaseCLRT_zBuffer(zbuffers[k]);
	keys.clear();

//  flushed_printf("kernels...\n");
	for (auto &p : kernels) keys.push_back(p.first);
	for (auto &k : keys) releaseCLRT_Kernel(kernels[k]);
	keys.clear();

//  flushed_printf("programs...\n");
	for (auto &p : programs) keys.push_back(p.first);
	for (auto &k : keys) releaseCLRT_Program(programs[k]);
	keys.clear();

//  flushed_printf("exeQueues...\n");
	for (auto &p : exeQueues) keys.push_back(p.first);
	for (auto &k : keys) releaseCLRT_Queue(exeQueues[k]);
	keys.clear();

//  flushed_printf("memory_queue_...\n");
	clReleaseCommandQueue(memory_queue_);
//  flushed_printf("context_...\n");
	clReleaseContext(context_);
//  flushed_printf("device_id_...\n");
	clReleaseDevice(device_id_);
//  flushed_printf("close opencl so...\n");
	close_opencl_so(&oclLibHandle_);

	// clear clInfo
	for (int kI = 0; kI < clInfo.num_platforms; ++kI) {
		for (int kJ = 0; kJ < clInfo.platforms[kI].num_devices; ++kJ) {
			delete clInfo.platforms[kI].devices[kJ].maxWGSizes;
		}
		delete clInfo.platforms[kI].devices;
	}
	delete clInfo.platforms;
}

CLRT_Queue CLRT::createSubDeviceExecutionProfilingQueue(uint Log2NumCUs, uint SDIndex) {
	auto queue = new InternalCmdQ_;
	queue->command_queue =
			clCreateCommandQueue(context_,
													 sub_deviceCollection[Log2NumCUs]->at(SDIndex),
													 CL_QUEUE_PROFILING_ENABLE,
													 &err,
													 oclLibHandle_);
	CLRT_ERR(err);
	exeQueues[queue->getTrackNum()] = queue;
	return queue;
}

CLRT_Queue CLRT::createCLExecutionQueue(cl_device_id Device, cl_command_queue_properties Properties) {
	auto queue = new InternalCmdQ_;
	queue->command_queue = clCreateCommandQueue(context_, Device, Properties, &err, oclLibHandle_);
	CLRT_ERR(err);
	exeQueues[queue->getTrackNum()] = queue;
	return queue;
}

CLRT_Queue CLRT::createCLExecutionQueue(cl_command_queue_properties Properties) {
	auto queue = new InternalCmdQ_;
	queue->command_queue = clCreateCommandQueue(context_, device_id_, Properties, &err, oclLibHandle_);
	CLRT_ERR(err);
	exeQueues[queue->getTrackNum()] = queue;
	return queue;
}

CLRT_Queue CLRT::createCLExecutionProfilingQueue() {
	return createCLExecutionQueue(CL_QUEUE_PROFILING_ENABLE);
}

void CLRT::pushBuffer(cl_mem CL_mem, size_t Offset, size_t nByte, void *HostPtr) {
	CLRT_ERR(clEnqueueWriteBuffer(memory_queue_, CL_mem, CL_TRUE,
																Offset, nByte, HostPtr,
																0, nullptr, nullptr, oclLibHandle_));
	clFinish(memory_queue_, oclLibHandle_);
}

void CLRT::pullBuffer(cl_mem CL_mem, size_t Offset, size_t nByte, void *HostPtr) {
	CLRT_ERR(clEnqueueReadBuffer(memory_queue_, CL_mem, CL_TRUE,
															 Offset, nByte, HostPtr,
															 0, nullptr, nullptr, oclLibHandle_));
	clFinish(memory_queue_, oclLibHandle_);
}

CLRT_Buffer CLRT::createBuffer(cl_mem_flags Flags, size_t nByte, void *HostPtr) {
	auto buffer = new InternalBuffer_(nByte);
	buffer->SIGNATURE = BUF;
	buffer->sizeByte = nByte;
	if (HostPtr)
		buffer->mem = clCreateBuffer(context_, Flags, nByte, HostPtr, &err, oclLibHandle_);
	else
		buffer->mem = clCreateBuffer(context_, Flags, nByte, nullptr, &err, oclLibHandle_);
	CLRT_ERR(err);
	buffers[buffer->getTrackNum()] = buffer;
	return buffer;
}

uint CLRT::getMaxSubDevice(int PlatformIdx, int DeviceIdx) const {
	return this->clInfo.platforms[PlatformIdx].devices[DeviceIdx].maxSubDevices;
}

void CLRT::
initSubDevices() {
	uint numCUPerDevice = getMaxSubDevice(0, 0);
	if (numCUPerDevice <= 1) { // device doesn't support splitting sub-devices
		this->numSubDeviceDivisions = 0;
		return;
	}
	while (numCUPerDevice > 0) {
		numSubDeviceDivisions++;
		numCUPerDevice = numCUPerDevice >> 1;
	}
	numCUPerDevice = 1; // reset to 1
	uint subDeviceDivisions[numSubDeviceDivisions];
//  uint *subDeviceDivisions = (uint *)calloc(numSubDeviceDivisions, sizeof(uint));
	for (int kI = 0; kI < numSubDeviceDivisions; ++kI) {
		subDeviceDivisions[kI] = numCUPerDevice << kI; // each sub-device can have 2^kI CUs
	}
	cl_device_partition_property props[3];
	props[0] = CL_DEVICE_PARTITION_EQUALLY;  // Equally
	props[2] = 0;                            // End of the property list

	for (int i = 0; i < numSubDeviceDivisions; ++i) {
		numCUPerDevice = subDeviceDivisions[i];
		props[1] = numCUPerDevice;     // compute units per sub-device
		// Specifies the size of the out_devices array: max num of sub-devices(num cpus) / num of cu(cpu core) per device
		uint numSubDevices = getMaxSubDevice(0, 0) / numCUPerDevice;
// Provides a buffer for the generated subdevices with a number
// of elements specified by num_sub_devices
		sub_deviceCollection.push_back(new std::vector<cl_device_id>(numSubDevices));
// clCreateSubDevices returns the number of subdevices
// in which the device may be partitioned into considering the
// partition type and the other values specified in the property list
		uint numSubDeviceRet = 0;
		CLRT_ERR(clCreateSubDevices(device_id_,
																props,
																numSubDevices,
																sub_deviceCollection[i]->data(),
																&numSubDeviceRet, oclLibHandle_));
		uint minSubDev = MIN(numSubDevices, numSubDeviceRet);
		for (int j = 0; j < minSubDev; ++j) {
			cl_command_queue queue =
					clCreateCommandQueue(context_,
															 sub_deviceCollection[i]->at(j),
															 CL_QUEUE_PROFILING_ENABLE,
															 &err,
															 oclLibHandle_);
			CLRT_ERR(err);
		}
	}
}

cl_event CLRT::
createUserEvent() {
	cl_event event = clCreateUserEvent(context_, &err, oclLibHandle_);
	CLRT_ERR(err);
	return event;
}

void CLRT::
setUserEventComplete(cl_event &event) const {
	if (event == nullptr) {
		std::cerr << "Setting a null event complete!\n";
		exit(-1);
	}
	CLRT_ERR(clSetUserEventStatus(event, CL_COMPLETE, oclLibHandle_));
}

void CLRT::
releaseEvent(cl_event &event) const {
	if (event == nullptr) return;
	CLRT_ERR(clReleaseEvent(event, oclLibHandle_));
	event = nullptr;
}

void CLRT::releaseCLRT_Queue(CLRT_Queue queue) {
	if (exeQueues.count(queue->getTrackNum()) > 0) {
		clReleaseCommandQueue(queue->command_queue);
		exeQueues.erase(queue->getTrackNum());
		delete queue;
	}
}
