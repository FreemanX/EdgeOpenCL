#include <unistd.h>
#include "EDCL.h"
#include "testlib_cpu_vector_add.h"
#include "testlib_macros.h"
#include "utils.h"

#define FLOAT_PTR(EDBuffer) HOST_PTR(float, EDBuffer)

void vecAdd_test(EDCL *edcl);

int main() {
  EDCL *edcl = new EDCL();
  vecAdd_test(edcl);
  delete edcl;
  return 0;
}

void vecAdd_test(EDCL *edcl) {
  const char *kernelSource =
#include "vecAdd.cl"

//#define PRINT_MAT
#ifdef PRINT_MAT
  unsigned int vDim = 4;
#else
  unsigned int vDim = 4096;
#endif
  uint size = vDim * vDim;
  size_t nByte = size * sizeof(float);
  double timer;
  int cmpResult;

  // Init buffers
  EDBuffer bufferA = edcl->createBuffer(nByte);
  EDBuffer bufferB = edcl->createBuffer(nByte);
  TIME_CODE_ONCE(
	  timer,
	  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferA), size);
		  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferB), size);
  )
  EDBuffer bufferC = edcl->createBuffer(nByte);
  auto *c_verify = static_cast<float *>(malloc(nByte));
  memset(bufferC->hostPtr, 0, nByte);
  memset(c_verify, 0, nByte);
  vecAddOptimized(FLOAT_PTR(bufferA), FLOAT_PTR(bufferB), c_verify, size);
#ifdef PRINT_MAT
  flushed_printf("A:\n");
  printMatrix(FLOAT_PTR(vecA), vDim, vDim);
  flushed_printf("B:\n");
  printMatrix(FLOAT_PTR(bufferB), vDim, vDim);
  flushed_printf("ground truth C:\n");
  printMatrix(c_verify, vDim, vDim);
#endif

  // Init kernels
  EDProgram program = edcl->createProgram(&kernelSource);
  EDKernel kernel = edcl->createKernel(program, "vecAdd", 4);
  size_t gw = size;
  kernel->configKernel(1, &gw, nullptr, bufferA, bufferB, bufferC, &size);

// Init execution
  EDQueue cpuQ = edcl->createDeviceCmdQueueProfilingEnabled(CPUQueue);
  edcl->confirmExeEnv(cpuQ, kernel);

  // Execute
  edcl->executeKernel(cpuQ, kernel, 0, nullptr, nullptr);
  cpuQ->finish();

  cmpResult = memcmp(bufferC->hostPtr, c_verify, nByte);
  flushed_printf("CPU vecAdd done, memcmp: %d, %s\n", cmpResult, cmpResult == 0 ? "PASS" : "FAIL");
#ifdef PRINT_MAT
  flushed_printf("CPU out C(%p): \n", bufferC->hostPtr);
  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
#endif
  memset(bufferC->hostPtr, 0, nByte);
  for (int kI = 0; kI < size; ++kI) FLOAT_PTR(bufferC)[kI] = (float)kI;
#ifdef PRINT_MAT
  flushed_printf("After clear bufferC(%p): \n", bufferC->hostPtr);
  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
#endif

  EDQueue gpuQ = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  edcl->confirmExeEnvAndExecuteSingleKernel(gpuQ, kernel, 0, nullptr, nullptr);
  gpuQ->finish();
  cmpResult = memcmp(bufferC->hostPtr, c_verify, nByte);
  flushed_printf("GPU vecAdd done, memcmp: %d, %s\n", cmpResult, cmpResult == 0 ? "PASS" : "FAIL");
#ifdef PRINT_MAT
  flushed_printf("GPU out C(%p): \n", bufferC->hostPtr);
  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
#endif

//  edcl->releaseEDBuffer(bufferA);
//  edcl->releaseEDBuffer(bufferB);
//  edcl->releaseEDBuffer(bufferC);
  edcl->releaseEDQueue(cpuQ);
  edcl->releaseEDQueue(gpuQ);
  delete kernel;
  edcl->releaseEDProgram(program);
}

