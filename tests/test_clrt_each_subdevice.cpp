#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include "CLRT.h"
#include "utils.h"
#include <vector>

const char *source =
#include "vecAdd.cl"

void testOpenCLLibrary(const char *path);

void vecAddOptimized(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddNoneOptimized(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP8Opt(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP4Opt(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP8(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP4(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP1(const float *va, const float *vb, float *vc, unsigned int n);

#define LOOP_NUM 50

#define TIME_CODE_LOOP(TIMER, CODE_BLOCK) \
{ \
    TIMER=getCurrentTime(); \
    for (int loopCnt = 0; loopCnt < LOOP_NUM; ++loopCnt) \
        {CODE_BLOCK} \
    TIMER=(getCurrentTime()-TIMER)/LOOP_NUM; \
}

#define TIME_CODE_ONCE(TIMER, CODE_BLOCK) \
{ \
    TIMER=getCurrentTime(); \
        {CODE_BLOCK} \
    TIMER=(getCurrentTime()-TIMER); \
}

float *a;
float *b;
float *c;
float *c_cl;
unsigned int vDim = 4096;
unsigned int size = vDim * vDim;
size_t nByte = size * sizeof(float);
double timer;

int main() {
  flushed_printf("Initializing vector(%d)... ", size);
  TIME_CODE_ONCE(
	  timer,
	  a = (float *)malloc(nByte);
		  b = (float *)malloc(nByte);
		  randArrayGenerator<float>(0.001, 1, a, size);
		  randArrayGenerator<float>(0.001, 1, b, size);
		  c = (float *)malloc(nByte);
		  c_cl = (float *)malloc(nByte);
  )
  flushed_printf("done. Time cost: %lf\n", timer);

  flushed_printf("Native CPU: ");
  TIME_CODE_ONCE(timer, vecAddNoneOptimized(a, b, c, size);)
  flushed_printf("\t\t%lf sec. \n", timer);

  flushed_printf("Native CPU Optimized: ");
  TIME_CODE_ONCE(timer, vecAddOptimized(a, b, c, size);)
  flushed_printf("\t%lf sec. \n", timer);

  flushed_printf("1 threads omp: ");
  TIME_CODE_ONCE(timer, vecAddOMP1(a, b, c, size);)
  flushed_printf("\t\t%lf sec. \n", timer);

  flushed_printf("4 threads omp: ");
  TIME_CODE_ONCE(timer, vecAddOMP4(a, b, c, size);)
  flushed_printf("\t\t%lf sec. \n", timer);

  flushed_printf("8 threads omp: ");
  TIME_CODE_ONCE(timer, vecAddOMP8(a, b, c, size);)
  flushed_printf("\t\t%lf sec. \n", timer);

  flushed_printf("4 threads omp_opt: ");
  TIME_CODE_ONCE(timer, vecAddOMP4Opt(a, b, c, size);)
  flushed_printf("\t%lf sec. \n", timer);

  flushed_printf("8 threads omp_opt: ");
  TIME_CODE_ONCE(timer, vecAddOMP8Opt(a, b, c, size);)
  flushed_printf("\t%lf sec. \n", timer);

  flushed_printf("CPU: \n");
  testOpenCLLibrary(DEFAULT_POCL_PATH);
  flushed_printf("GPU: \n");
  testOpenCLLibrary(DEFAULT_GPU_OCL_PATH);
  return 0;
}

void exeVecAdd(CLRT *clrt, CLRT_Program libVecAdd, cl_command_queue exeQ) {
  cl_event event;
  // Test with InternalBuffer_
  auto bufferA = clrt->createBuffer(nByte, a);
  auto bufferB = clrt->createBuffer(nByte, nullptr);
  clrt->pushBuffer(bufferB, 0, nByte, b);
  auto bufferC = clrt->createBuffer(nByte, nullptr);
  auto vecAdd = clrt->createKernel(libVecAdd, "vecAdd", 4); // Kernel -> function in library
  CLRT_SET_ARGS(vecAdd, bufferA, bufferB, bufferC, &size)
  size_t globalSize = size;

  printf("  OpenCL vecAdd:");
  //Dry run
  clrt->execKernel(exeQ, vecAdd, 1, nullptr, &globalSize, nullptr, 0, nullptr, nullptr);
  clFinish(exeQ);
  TIME_CODE_LOOP (
	  timer,
	  clrt->execKernel(exeQ, vecAdd, 1, nullptr, &globalSize, nullptr, 0, nullptr, &event);
		  clWaitForEvents(1, &event, clrt->oclLibHandle_);
  )
  clFinish(exeQ);
  clrt->pullBuffer(bufferC, 0, size * sizeof(float), c_cl);
//  flushed_printf(" %lf sec. OpenCL event timer: %lf sec, ",
//				 timer, CLRT_GET_EVENT_TIME(&event));
  flushed_printf("\t%lf sec. ", timer);
  flushed_printf(" memcmp = %d \n", memcmp(c, c_cl, nByte));

//  flushed_printf("Cleaning... ");
//  flushed_printf("vecAdd ");
  clrt->releaseCLRT_Kernel(vecAdd);
//  flushed_printf("bufferA ");
  clrt->releaseCLRT_Buffer(bufferA);
//  flushed_printf("bufferB ");
  clrt->releaseCLRT_Buffer(bufferB);
//  flushed_printf("bufferC ");
  clrt->releaseCLRT_Buffer(bufferC);
  memset(c_cl, 0, nByte);
//  flushed_printf("Done! \n");
}

void exeMatAdd(CLRT *clrt, CLRT_Program libVecAdd, cl_command_queue exeQ) {
  cl_event event;
  // Test with zBuffer
  auto matAdd = clrt->createKernel(libVecAdd, "matAdd", 5);
  auto matAdd_ArgTester = clrt->createKernel(libVecAdd, "matAdd", 5);
  auto zbufferA = clrt->createZBuffer(nByte, a);
  auto zbufferB = clrt->createZBuffer(nByte, nullptr);
  auto zbufferC = clrt->createZBuffer(nByte, nullptr);
  clrt->mapZBuffer(zbufferB);
  clrt->mapZBuffer(zbufferC);
  memcpy(zbufferA->hostPtr, a, nByte);
  memcpy(zbufferB->hostPtr, b, nByte);
  memset(zbufferC->hostPtr, 0, nByte);
  clrt->unmapZBuffer(zbufferA);
  clrt->unmapZBuffer(zbufferB);
  clrt->unmapZBuffer(zbufferC);
  CLRT_SET_ARGS(matAdd_ArgTester, zbufferA, zbufferB, zbufferC, &vDim, &vDim)
  clrt->setKernelArgWithLists(matAdd->kernel_,
							  matAdd_ArgTester->argTypeList,
							  matAdd_ArgTester->argPtrList,
							  matAdd_ArgTester->argSizeList,
							  matAdd->numArgs);
  clrt->releaseCLRT_Kernel(matAdd_ArgTester);
  size_t gw[] = {vDim / 4, vDim};
  printf("  OpenCL matAdd:");
  clrt->execKernel(exeQ, matAdd, 2, nullptr, gw, nullptr, 0, nullptr, nullptr);
  clFinish(exeQ);
  TIME_CODE_LOOP(
	  timer,
	  clrt->execKernel(exeQ, matAdd, 2, nullptr, gw, nullptr, 0, nullptr, &event);
		  clWaitForEvents(1, &event, clrt->oclLibHandle_);
  )
  clFinish(exeQ);
  clrt->mapZBuffer(zbufferC);
  flushed_printf("\t%lf sec. ", timer);
//  flushed_printf(" %lf sec. OpenCL event timer: %lf sec, ",
//				 timer, CLRT_GET_EVENT_TIME(&event));
  flushed_printf("memcmp = %d \n", memcmp(c, zbufferC->hostPtr, nByte));

  // Test release functions
//  flushed_printf("Cleaning... ");
//  flushed_printf("matAdd ");
  clrt->releaseCLRT_Kernel(matAdd);
//  flushed_printf("zbufferA ");
  clrt->releaseCLRT_zBuffer(zbufferA);
//  flushed_printf("zbufferB ");
  clrt->releaseCLRT_zBuffer(zbufferB);
//  flushed_printf("zbufferC ");
  clrt->releaseCLRT_zBuffer(zbufferC);
}

void testOpenCLLibrary(const char *path) {
  CLRT *clrt = new CLRT(path);
  auto libVecAdd = clrt->createProgramFromSource(&source); // program -> library
  cl_command_queue exeQ;
  if (strcmp(path, DEFAULT_POCL_PATH) == 0) {
	uint cuPerDevice[4] = {1, 2, 4, 8};
	cl_device_partition_property props[3];
	props[0] = CL_DEVICE_PARTITION_EQUALLY;  // Equally
	props[2] = 0;                            // End of the property list

	for (unsigned int kI : cuPerDevice) {
	  props[1] = kI;                            // compute units per sub-device
	  printf(" Subdevices, cu per device: %d\n", kI);
	  // Specifies the size of the out_devices array
	  cl_uint num_sub_devices = 8 / kI;
// Provides a buffer for the generated subdevices with a number
// of elements specified by num_sub_devices
	  cl_device_id sub_device_ids[num_sub_devices];
// clCreateSubDevices returns the number of subdevices
// in which the device may be partitioned into considering the
// partition type and the other values specified in the property list
	  cl_uint num_devices_ret = 0;

// Create the sub-devices for the device_id device
	  cl_int err;
	  clCreateSubDevices(clrt->device_id_, props, num_sub_devices, sub_device_ids, &num_devices_ret);
	  for (int kJ = 0; kJ < num_sub_devices; ++kJ) {
		flushed_printf("Running on SD[%d]: ", kJ);
		exeQ = clCreateCommandQueue(clrt->context_, sub_device_ids[kJ], CL_QUEUE_PROFILING_ENABLE, &err);
		CLRT_ERR(err);
		exeVecAdd(clrt, libVecAdd, exeQ);
//  flushed_printf("exeQ ");
		clReleaseCommandQueue(exeQ);
	  }
	}

  } else {
	exeQ = clrt->createCLExecutionProfilingQueue()->command_queue;
	exeVecAdd(clrt, libVecAdd, exeQ);
//  flushed_printf("exeQ ");
	clReleaseCommandQueue(exeQ);
  }

  if (strcmp(path, DEFAULT_POCL_PATH) == 0) {
	uint cuPerDevice[4] = {1, 2, 4, 8};
	cl_device_partition_property props[3];
	props[0] = CL_DEVICE_PARTITION_EQUALLY;  // Equally
	props[2] = 0;                            // End of the property list

	for (unsigned int kI : cuPerDevice) {
	  props[1] = kI;                            // compute units per sub-device
	  printf(" Subdevices, cu per device: %d\n", kI);
	  // Specifies the size of the out_devices array
	  cl_uint num_sub_devices = 8 / kI;
// Provides a buffer for the generated subdevices with a number
// of elements specified by num_sub_devices
	  cl_device_id sub_device_ids[num_sub_devices];
// clCreateSubDevices returns the number of subdevices
// in which the device may be partitioned into considering the
// partition type and the other values specified in the property list
	  cl_uint num_devices_ret = 0;

// Create the sub-devices for the device_id device
	  cl_int err;
	  clCreateSubDevices(clrt->device_id_, props, num_sub_devices, sub_device_ids, &num_devices_ret);
	  for (int kJ = 0; kJ < num_sub_devices; ++kJ) {
		flushed_printf("Running on SD[%d]: ", kJ);
		exeQ = clCreateCommandQueue(clrt->context_, sub_device_ids[kJ], CL_QUEUE_PROFILING_ENABLE, &err);
		CLRT_ERR(err);
		exeMatAdd(clrt, libVecAdd, exeQ);
//  flushed_printf("exeQ ");
		clReleaseCommandQueue(exeQ);
	  }
	}

  } else {
	exeQ = clrt->createCLExecutionProfilingQueue()->command_queue;
	exeMatAdd(clrt, libVecAdd, exeQ);
//  flushed_printf("exeQ ");
	clReleaseCommandQueue(exeQ);
  }


//  flushed_printf("libVecAdd ");
  clrt->releaseCLRT_Program(libVecAdd);
//  flushed_printf("&clrt ");
  delete clrt;
//  flushed_printf("Done! \n");
}

void vecAddNoneOptimized(const float *va,
						 const float *vb,
						 float *vc,
						 unsigned int n)__attribute__((optnone)) {
  int i;
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOptimized(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP8Opt(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(8)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP4Opt(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(4)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP8(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(8)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP4(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(4)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP1(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(1)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}
