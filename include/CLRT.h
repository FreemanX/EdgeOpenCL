// CL runtime
#pragma ide diagnostic ignored "OCUnusedMacroInspection"
#ifndef EDCL_INCLUDE_CLRT_H_
#define EDCL_INCLUDE_CLRT_H_
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "CL/cl.h"
#include "BufferDef.h"
#include "custom_cl.h"

typedef struct InternalKernel_ : public Trackable {
  void *oclLibHandle;     // internal
  std::string kernelName; // external
  cl_kernel kernel_;      // internal
  unsigned int numArgs;      // external
  unsigned int argIndex = 0; // internal used when auto set arguments

  /* argTypeList is used for determine weather to use sizeof(cl_mem) or
   * argSize when calling clSetKernelArg().
   * Pointers in argPtrList can be passed to clSetKernelArg() directly.
   * argSizeList contains the actual data size if arg is cl_mem. Used to
   * tell others how much memory is needed. */
  // argTypeList: one of BufferType, 2 categories: cl_mem or not cl_mem
  BufferType *argTypeList; // external
  // argPtrList: pointers point to cl_mem or a primitive type
  void **argPtrList; // external
  // argSizeList: the actual data size, not sizeof(cl_mem) if arg is cl_mem
  size_t *argSizeList; // external
} InternalKernel_;

typedef struct InternalProgram_ : public Trackable {
  cl_program program_{}; // internal
} InternalProgram_;

typedef struct InternalCmdQ_ : public Trackable {
  cl_command_queue command_queue{}; // internal
} InternalCmdQ_;

typedef struct {
  cl_device_id deviceId;
  int num_CUs;
  int localMemorySizeB;
  int WGDimSize;
  size_t *maxWGSizes;
  int maxAllocSizeB;
  uint maxSubDevices;
  // more device info adds here...
} InternalDevice_CL_;

typedef struct {
  cl_platform_id platformId;
  cl_uint num_devices;
  InternalDevice_CL_ *devices;
} InternalPlatform_CL_;

typedef struct {
  cl_uint num_platforms;
  InternalPlatform_CL_ *platforms;
} InternalInfo_CL_;

/*-------------------------*/
/* macros, types for users */
/*-------------------------*/
typedef struct InternalProgram_ *CLRT_Program;
typedef struct InternalKernel_ *CLRT_Kernel;
typedef struct InternalCmdQ_ * CLRT_Queue;
typedef struct InternalBuffer_ *CLRT_Buffer;   // normal CL buffer: cl_mem
typedef struct InternalZBuffer_ *CLRT_zBuffer; // zero copy buffer: cl_mem + host_ptr



#define DEFAULT_POCL_PATH "/data/local/tmp/pocl/lib/libOpenCL.so.2.6.0"
#define DEFAULT_GPU_OCL_PATH "/vendor/lib64/libOpenCL.so"
#define CLRT_ERR(x) CLRT::CLRTChkErr(x)

// Load OpenCL strings like #include "<source file>.cl"
// Reference: https://stackoverflow.com/a/5034490/12709772
// https://github.com/FreemanX/clpeak/blob/master/src/clpeak.cpp#L6-L11
// Note that Stringifying requires a new line after hash defines:
// 		https://github.com/FreemanX/clpeak/blob/master/src/kernels/compute_dp_kernels.cl
#define MSTRINGIFY(...) #__VA_ARGS__
/* user macros ------------ */


// OpenCL Runtime
class CLRT {

 public: // properties
  cl_platform_id platform_id_{};
  cl_device_id device_id_{};
  cl_context context_{};
  void *oclLibHandle_{};
  InternalInfo_CL_ clInfo{};
  // sub-devices and corresponding queues
  uint numSubDeviceDivisions = 0; // How many sub-device options?
// Array of sub-device vectors. Each vector contains sub-devices with the same number of CUs.
  std::vector<std::vector<cl_device_id> *> sub_deviceCollection;

 private: // properties
  std::string oclLibPath_;
  cl_command_queue memory_queue_{};
  cl_int err{};
  std::unordered_map<int, CLRT_Queue> exeQueues;
  std::unordered_map<int, CLRT_Program> programs;
  std::unordered_map<int, CLRT_Kernel> kernels;
  std::unordered_map<int, CLRT_Buffer> buffers;
  std::unordered_map<int, CLRT_zBuffer> zbuffers;

 public: // functions
  std::string getOCLLibPath() { return this->oclLibPath_; }

  CLRT();

  explicit CLRT(std::string OclLibPath);

  static CLRT *createCPURunTime() { return new CLRT(DEFAULT_POCL_PATH); }

  static CLRT *createGPURunTime() { return new CLRT(DEFAULT_GPU_OCL_PATH); }

  CLRT_Queue createCLExecutionQueue(cl_device_id Device, cl_command_queue_properties Properties);

  CLRT_Queue createCLExecutionQueue(cl_command_queue_properties Properties);

  CLRT_Queue createCLExecutionProfilingQueue();

  CLRT_Queue createSubDeviceExecutionProfilingQueue(uint Log2NumCUs, uint SDIndex);

  static std::string CLRTErrMsg(cl_int error);

  static void CLRTChkErr(cl_int error);

  void printCLInfo() const;

  uint getMaxSubDevice(int PlatformIdx, int DeviceIdx) const;

  double getExeTime(cl_event const *e) const;

  CLRT_Program createProgramFromSource(const char **ProgramSource);

  CLRT_Program createProgramFromSourceWithOptions(const char **ProgramSource, const char *options);

  CLRT_Kernel createKernel(CLRT_Program &Program,
						   std::string Name,
						   unsigned int NumArgs);

  // Set to static in order to use macro
  static void setKernelArgUsedByMacro_CLRT(CLRT_Kernel Kernel,
										   size_t ArgObjSize,
										   void *ArgPtr);

  void setKernelArgAtIndex(CLRT_Kernel Kernel,
						   uint Index,
						   BufferType ArgType,
						   size_t ArgSize,
						   void *ArgPtr) const;

  void setKernelArgWithLists(cl_kernel ClKernel,
							 const BufferType *TypeList,
							 void **PtrList,
							 size_t *SizeList,
							 uint NumArgs) const;

  void setKernelArgWithLists(CLRT_Kernel CLRTKernel,
							 const BufferType *TypeList,
							 void **PtrList,
							 size_t *SizeList,
							 uint NumArgs) const;

  void execKernel(cl_command_queue ExeQ,
				  CLRT_Kernel Kernel,
				  unsigned int WorkDim,
				  const size_t *Offset,
				  const size_t *GlobalSize,
				  const size_t *LocalSize,
				  unsigned int NumEventsWait,
				  const cl_event *EventWaitList,
				  cl_event *event) const;

  CLRT_Buffer createBuffer(size_t nByte, void *HostPtr);

  CLRT_Buffer createBuffer(cl_mem_flags Flags, size_t nByte, void *HostPtr);

  void pushBuffer(CLRT_Buffer Buffer, size_t Offset, size_t nByte, void *HostPtr);

  void pullBuffer(CLRT_Buffer Buffer, size_t Offset, size_t nByte, void *HostPtr);

  void pushBuffer(cl_mem CL_mem, size_t Offset, size_t nByte, void *HostPtr);

  void pullBuffer(cl_mem CL_mem, size_t Offset, size_t nByte, void *HostPtr);

  // HostPtr: the pointer to the buffer data that may already be allocated by the application.
  // zBuffer's host_ptr is mapped by default
  CLRT_zBuffer createZBuffer(size_t nByte, void *HostPtr);

  void *mapZBuffer(CLRT_zBuffer ZBuffer);

  void unmapZBuffer(CLRT_zBuffer ZBuffer);

  cl_event createUserEvent();

  void setUserEventComplete(cl_event &event) const;

  void releaseEvent(cl_event &event) const;

  void releaseCLRT_Program(CLRT_Program Program);

  void releaseCLRT_Kernel(CLRT_Kernel Kernel);

  void releaseCLRT_Buffer(CLRT_Buffer Buffer);

  void releaseCLRT_zBuffer(CLRT_zBuffer zBuffer);

  void releaseCLRT_Queue(CLRT_Queue queue);

  /* User can choose to release by their own
   * or everything will be cleared when
   * clrt object is deleted*/
  ~CLRT();

  // public functions end
 private: // functions
  void initRT(cl_platform_id Platform, cl_device_id Device);
  void queryCLInfo(InternalInfo_CL_ *CLInfo) const;
  cl_program createProgramFromSource_cl(cl_context Context,
										cl_device_id *Devices,
										unsigned int NumDevices,
										const char **Source,
										const char *Options);

  void initSubDevices();
  // private functions end
  cl_kernel createKernel_cl(cl_program Program, std::string &Name);
};

// ------------
// CL Kernel args handle
// Reference: https://github.com/champyen/cltk/blob/master/include/cltk.h
// ------------
#define CLRT_ARGS_HANDLE(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _A, _B, _C, _D, \
                     _E, _F, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, \
                     _1A, _1B, _1C, _1D, _1E, _1F, N, ...)                     \
  N

#define CLRT_SET_ARG(Kernel, ArgPtr)                                                  \
  { CLRT::setKernelArgUsedByMacro_CLRT(Kernel, sizeof(*ArgPtr), ArgPtr); }

#define _ARG_0(...)
#define _ARG_1(Kernel, ArgPtr)      CLRT_SET_ARG(Kernel, ArgPtr)
#define _ARG_2(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1( Kernel, __VA_ARGS__)
#define _ARG_3(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_2( Kernel, __VA_ARGS__)
#define _ARG_4(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_3( Kernel, __VA_ARGS__)
#define _ARG_5(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_4( Kernel, __VA_ARGS__)
#define _ARG_6(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_5( Kernel, __VA_ARGS__)
#define _ARG_7(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_6( Kernel, __VA_ARGS__)
#define _ARG_8(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_7( Kernel, __VA_ARGS__)
#define _ARG_9(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_8( Kernel, __VA_ARGS__)
#define _ARG_A(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_9( Kernel, __VA_ARGS__)
#define _ARG_B(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_A( Kernel, __VA_ARGS__)
#define _ARG_C(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_B( Kernel, __VA_ARGS__)
#define _ARG_D(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_C( Kernel, __VA_ARGS__)
#define _ARG_E(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_D( Kernel, __VA_ARGS__)
#define _ARG_F(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_E( Kernel, __VA_ARGS__)
#define _ARG_10(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_F( Kernel, __VA_ARGS__)
#define _ARG_11(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_10(Kernel, __VA_ARGS__)
#define _ARG_12(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_11(Kernel, __VA_ARGS__)
#define _ARG_13(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_12(Kernel, __VA_ARGS__)
#define _ARG_14(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_13(Kernel, __VA_ARGS__)
#define _ARG_15(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_14(Kernel, __VA_ARGS__)
#define _ARG_16(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_15(Kernel, __VA_ARGS__)
#define _ARG_17(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_16(Kernel, __VA_ARGS__)
#define _ARG_18(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_17(Kernel, __VA_ARGS__)
#define _ARG_19(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_18(Kernel, __VA_ARGS__)
#define _ARG_1A(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_19(Kernel, __VA_ARGS__)
#define _ARG_1B(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1A(Kernel, __VA_ARGS__)
#define _ARG_1C(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1B(Kernel, __VA_ARGS__)
#define _ARG_1D(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1C(Kernel, __VA_ARGS__)
#define _ARG_1E(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1D(Kernel, __VA_ARGS__)
#define _ARG_1F(Kernel, ArgPtr, ...) CLRT_SET_ARG(Kernel, ArgPtr) _ARG_1E(Kernel, __VA_ARGS__)

// ... must all be pointers
#define CLRT_SET_ARGS(Kernel, ...)                                                     \
  {                                                                            \
    CLRT_ARGS_HANDLE("ignored", ##__VA_ARGS__, _ARG_1F, _ARG_1E, _ARG_1D, _ARG_1C, \
                 _ARG_1B, _ARG_1A, _ARG_19, _ARG_18, _ARG_17, _ARG_16,         \
                 _ARG_15, _ARG_14, _ARG_13, _ARG_12, _ARG_11, _ARG_10, _ARG_F, \
                 _ARG_E, _ARG_D, _ARG_C, _ARG_B, _ARG_A, _ARG_9, _ARG_8,       \
                 _ARG_7, _ARG_6, _ARG_5, _ARG_4, _ARG_3, _ARG_2, _ARG_1,       \
                 _ARG_0)                                                       \
    (Kernel, ##__VA_ARGS__)                                                    \
  }

#endif //EDCL_INCLUDE_CLRT_H_
