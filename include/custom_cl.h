#ifndef EDCL_INCLUDE_CUSTOM_CL_H_
#define EDCL_INCLUDE_CUSTOM_CL_H_
#include "CL/cl.h"
#include "libopencl.h"
#include <dlfcn.h>
#include <cassert>

static cl_int clWaitForEvents(cl_uint num_events, const cl_event *event_list, void *handle) {
  f_clWaitForEvents func;

  func = (f_clWaitForEvents)dlsym(handle, "clWaitForEvents");
  if (func) {
	return func(num_events, event_list);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clGetEventProfilingInfo(cl_event event,
									  cl_profiling_info param_name,
									  size_t param_value_size,
									  void *param_value,
									  size_t *param_value_size_ret,
									  void *so_handle) {
  f_clGetEventProfilingInfo func;

  func = (f_clGetEventProfilingInfo)dlsym(so_handle, "clGetEventProfilingInfo");
  if (func) {
	return func(event, param_name, param_value_size, param_value, param_value_size_ret);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clEnqueueNDRangeKernel(cl_command_queue command_queue,
									 cl_kernel kernel,
									 cl_uint work_dim,
									 const size_t *global_work_offset,
									 const size_t *global_work_size,
									 const size_t *local_work_size,
									 cl_uint num_events_in_wait_list,
									 const cl_event *event_wait_list,
									 cl_event *event,
									 void *so_handle) {
  f_clEnqueueNDRangeKernel func;

  func = (f_clEnqueueNDRangeKernel)dlsym(so_handle, "clEnqueueNDRangeKernel");
  if (func) {
	return func(command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size,
				num_events_in_wait_list, event_wait_list, event);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clReleaseEvent(cl_event event, void *so_handle) {
  f_clReleaseEvent func;

  func = (f_clReleaseEvent)dlsym(so_handle, "clReleaseEvent");
  if (func) {
	return func(event);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clSetUserEventStatus(cl_event event, cl_int execution_status, void *so_handle) {
  f_clSetUserEventStatus func;

  func = (f_clSetUserEventStatus)dlsym(so_handle, "clSetUserEventStatus");
  if (func) {
	return func(event, execution_status);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_event clCreateUserEvent(cl_context context, cl_int *errcode_ret, void *so_handle) {
  f_clCreateUserEvent func;

  func = (f_clCreateUserEvent)dlsym(so_handle, "clCreateUserEvent");
  if (func) {
	return func(context, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_command_queue clCreateCommandQueue(cl_context context,
											 cl_device_id device,
											 cl_command_queue_properties properties,
											 cl_int *errcode_ret, void *so_handle) {
  f_clCreateCommandQueue func;
  func = (f_clCreateCommandQueue)dlsym(so_handle, "clCreateCommandQueue");
  if (func) {
	return func(context, device, properties, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_int clCreateSubDevices(cl_device_id in_device,
								 const cl_device_partition_property *properties,
								 cl_uint num_devices,
								 cl_device_id *out_devices,
								 cl_uint *num_devices_ret, void *so_handle) {
  f_clCreateSubDevices func;
  func = (f_clCreateSubDevices)dlsym(so_handle, "clCreateSubDevices");
  if (func) {
	return func(in_device, properties, num_devices, out_devices, num_devices_ret);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_mem
clCreateBuffer(cl_context context,
			   cl_mem_flags flags,
			   size_t size,
			   void *host_ptr,
			   cl_int *errcode_ret,
			   void *so_handle) {
  f_clCreateBuffer func;

  func = (f_clCreateBuffer)dlsym(so_handle, "clCreateBuffer");
  if (func) {
	return func(context, flags, size, host_ptr, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_int clFinish(cl_command_queue command_queue, void *so_handle) {
  assert(so_handle != nullptr);
  f_clFinish func;

  func = (f_clFinish)dlsym(so_handle, "clFinish");
  if (func) {
	return func(command_queue);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clEnqueueReadBuffer(cl_command_queue command_queue,
								  cl_mem buffer,
								  cl_bool blocking_read,
								  size_t offset,
								  size_t size,
								  void *ptr,
								  cl_uint num_events_in_wait_list,
								  const cl_event *event_wait_list,
								  cl_event *event, void *so_handle) {
  f_clEnqueueReadBuffer func;

  func = (f_clEnqueueReadBuffer)dlsym(so_handle, "clEnqueueReadBuffer");
  if (func) {
	return func(command_queue, buffer, blocking_read, offset, size, ptr,
				num_events_in_wait_list, event_wait_list, event);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clEnqueueWriteBuffer(cl_command_queue command_queue,
								   cl_mem buffer,
								   cl_bool blocking_write,
								   size_t offset,
								   size_t size,
								   const void *ptr,
								   cl_uint num_events_in_wait_list,
								   const cl_event *event_wait_list,
								   cl_event *event, void *so_handle) {
  f_clEnqueueWriteBuffer func;

  func = (f_clEnqueueWriteBuffer)dlsym(so_handle, "clEnqueueWriteBuffer");
  if (func) {
	return func(command_queue, buffer, blocking_write, offset, size, ptr,
				num_events_in_wait_list, event_wait_list, event);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int
clReleaseMemObject(cl_mem memobj, void *so_handle) {
  assert(so_handle != nullptr);
  f_clReleaseMemObject func;
  func = (f_clReleaseMemObject)dlsym(so_handle, "clReleaseMemObject");
  if (func) {
	return func(memobj);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clEnqueueUnmapMemObject(cl_command_queue command_queue,
									  cl_mem memobj,
									  void *mapped_ptr,
									  cl_uint num_events_in_wait_list,
									  const cl_event *event_wait_list,
									  cl_event *event, void *so_handle) {
  f_clEnqueueUnmapMemObject func;

  func = (f_clEnqueueUnmapMemObject)dlsym(so_handle, "clEnqueueUnmapMemObject");
  if (func) {
	return func(command_queue, memobj, mapped_ptr, num_events_in_wait_list, event_wait_list, event);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static void *clEnqueueMapBuffer(cl_command_queue command_queue,
								cl_mem buffer,
								cl_bool blocking_map,
								cl_map_flags map_flags,
								size_t offset,
								size_t size,
								cl_uint num_events_in_wait_list,
								const cl_event *event_wait_list,
								cl_event *event,
								cl_int *errcode_ret,
								void *so_handle) {
  f_clEnqueueMapBuffer func;
  func = (f_clEnqueueMapBuffer)dlsym(so_handle, "clEnqueueMapBuffer");
  if (func) {
	return func(command_queue, buffer, blocking_map, map_flags, offset, size,
				num_events_in_wait_list, event_wait_list, event, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_int clSetKernelArg(cl_kernel kernel,
							 cl_uint arg_index,
							 size_t arg_size,
							 const void *arg_value, void *so_handle) {
  f_clSetKernelArg func;

  func = (f_clSetKernelArg)dlsym(so_handle, "clSetKernelArg");
  if (func) {
	return func(kernel, arg_index, arg_size, arg_value);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_program clCreateProgramWithSource(cl_context context,
											cl_uint count,
											const char **strings,
											const size_t *lengths,
											cl_int *errcode_ret, void *so_handle) {
  f_clCreateProgramWithSource func;

  func = (f_clCreateProgramWithSource)dlsym(so_handle, "clCreateProgramWithSource");
  if (func) {
	return func(context, count, strings, lengths, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_int clBuildProgram(cl_program program,
							 cl_uint num_devices,
							 const cl_device_id *device_list,
							 const char *options,
							 void (*pfn_notify)(__unused cl_program program, __unused void *user_data),
							 void *user_data, void *so_handle) {
  f_clBuildProgram func;
  func = (f_clBuildProgram)dlsym(so_handle, "clBuildProgram");
  if (func) {
	return func(program, num_devices, device_list, options, pfn_notify, user_data);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clGetProgramBuildInfo(cl_program program,
									cl_device_id device,
									cl_program_build_info param_name,
									size_t param_value_size,
									void *param_value,
									size_t *param_value_size_ret, void *so_handle) {
  f_clGetProgramBuildInfo func;

  func = (f_clGetProgramBuildInfo)dlsym(so_handle, "clGetProgramBuildInfo");
  if (func) {
	return func(program, device, param_name, param_value_size,
				param_value, param_value_size_ret);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_kernel clCreateKernel(cl_program program,
								const char *kernel_name,
								cl_int *errcode_ret, void *so_handle) {
  f_clCreateKernel func;

  func = (f_clCreateKernel)dlsym(so_handle, "clCreateKernel");
  if (func) {
	return func(program, kernel_name, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_context
clCreateContext(const cl_context_properties *properties,
				cl_uint num_devices,
				const cl_device_id *devices,
				void (*pfn_notify)(const char *, const void *, size_t, void *),
				void *user_data,
				cl_int *errcode_ret, void *so_handle) {
  f_clCreateContext func;

  func = (f_clCreateContext)dlsym(so_handle, "clCreateContext");
  if (func) {
	return func(properties, num_devices, devices, pfn_notify, user_data, errcode_ret);
  } else {
	return nullptr;
  }
}

static cl_int clGetPlatformIDs(cl_uint num_entries, cl_platform_id *platforms, cl_uint *num_platforms, void* so_handle) {
  f_clGetPlatformIDs func;

  func = (f_clGetPlatformIDs)dlsym(so_handle, "clGetPlatformIDs");
  if (func) {
	return func(num_entries, platforms, num_platforms);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clGetDeviceIDs(cl_platform_id platform,
					  cl_device_type device_type,
					  cl_uint num_entries,
					  cl_device_id *devices,
					  cl_uint *num_devices, void * so_handle) {
  f_clGetDeviceIDs func;

  func = (f_clGetDeviceIDs)dlsym(so_handle, "clGetDeviceIDs");
  if (func) {
	return func(platform, device_type, num_entries, devices, num_devices);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int
clGetDeviceInfo(cl_device_id device,
				cl_device_info param_name,
				size_t param_value_size,
				void *param_value,
				size_t *param_value_size_ret, void * so_handle) {
  f_clGetDeviceInfo func;

  func = (f_clGetDeviceInfo)dlsym(so_handle, "clGetDeviceInfo");
  if (func) {
	return func(device, param_name, param_value_size, param_value, param_value_size_ret);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int clGetPlatformInfo(cl_platform_id platform,
						 cl_platform_info param_name,
						 size_t param_value_size,
						 void *param_value,
						 size_t *param_value_size_ret, void * so_handle) {
  f_clGetPlatformInfo func;

  func = (f_clGetPlatformInfo)dlsym(so_handle, "clGetPlatformInfo");
  if (func) {
	return func(platform, param_name, param_value_size, param_value, param_value_size_ret);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int
clReleaseProgram(cl_program program, void *so_handle) {
  f_clReleaseProgram func;

  func = (f_clReleaseProgram)dlsym(so_handle, "clReleaseProgram");
  if (func) {
	return func(program);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int
clReleaseKernel(cl_kernel kernel, void * so_handle) {
  f_clReleaseKernel func;

  func = (f_clReleaseKernel)dlsym(so_handle, "clReleaseKernel");
  if (func) {
	return func(kernel);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

static cl_int
clReleaseDevice(cl_device_id device, void *so_handle) {
  f_clReleaseDevice func;

  func = (f_clReleaseDevice)dlsym(so_handle, "clReleaseDevice");
  if (func) {
	return func(device);
  } else {
	return CL_INVALID_PLATFORM;
  }
}

#endif //EDCL_INCLUDE_CUSTOM_CL_H_
