#include "EDCL.h"
#include <unistd.h>
#include "testlib_cpu_vector_add.h"
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

  unsigned int vDim = 1024;
  uint size = vDim * vDim;
  size_t nByte = size * sizeof(float);

  // Init buffers
  EDBuffer bufferA = edcl->createBuffer(nByte);
  EDBuffer bufferB = edcl->createBuffer(nByte);
  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferA), size);
  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferB), size);
  EDBuffer bufferC = edcl->createBuffer(nByte);
  auto *c_verify = static_cast<float *>(malloc(nByte));
  memset(bufferC->hostPtr, 0, nByte);
  memset(c_verify, 0, nByte);
  vecAddOptimized(FLOAT_PTR(bufferA), FLOAT_PTR(bufferB), c_verify, size);

  // Init kernels
  EDProgram program = edcl->createProgram(&kernelSource);
  EDKernel kernel = edcl->createKernel(program, "vecAdd", 4);
  size_t gw = size;
  kernel->configKernel(1, &gw, nullptr, bufferA, bufferB, bufferC, &size);

// Init execution
  EDQueue cpuQ = edcl->createDeviceCmdQueueProfilingEnabled(CPUQueue);
  edcl->confirmExeEnv(cpuQ, kernel);
//  cl_event event1;
//  edcl->executeKernel(cpuQ, kernel, 0, nullptr, &event1);

  cl_event event2;
  // Execute
  EDQueue gpuQ = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  edcl->confirmExeEnvAndExecuteSingleKernel(gpuQ, kernel, 0, nullptr, &event2);

  edcl->executeKernel(cpuQ, kernel, 1, &event2, nullptr);

  edcl->releaseEDBuffer(bufferA);
  edcl->releaseEDBuffer(bufferB);
  edcl->releaseEDBuffer(bufferC);
  edcl->releaseEDQueue(cpuQ);
  edcl->releaseEDQueue(gpuQ);
  delete kernel;
  edcl->releaseEDProgram(program);
}

