#include "EDCL.h"
#include "testlib_macros.h"

// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()

const char *source =
#include "vecAdd.cl"
uint vDim = 4096;
uint size = vDim * vDim;
size_t nByte = size * sizeof(float);
EDCL *edcl;
EDBuffer bufferA;
EDBuffer bufferB;
EDBuffer bufferC;
void splitKernelRun(EDKernel kernel, EDQueue exeQ, float);
__unused void printGlobalSize(EDKernel kernel);

int main() {
  edcl = new EDCL();

  bufferA = edcl->createBuffer(nByte);
  bufferB = edcl->createBuffer(nByte);
  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferA), size);
  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferB), size);
  bufferC = edcl->createBuffer(nByte);
  memset(bufferC->hostPtr, 0, nByte);

//  EDQueue exeQ = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  EDQueue exeQ = edcl->createSubDeviceExeQueue(2, 1);
  EDProgram program = edcl->createProgram(&source);

  EDKernel vecAdd = edcl->createKernel(program, "vecAdd", 4);
  size_t globalSize1D[] = {size, 1, 1};
  vecAdd->configKernel(1, globalSize1D, nullptr, bufferA, bufferB, bufferC, &size);

  EDKernel matAdd = edcl->createKernel(program, "matAdd", 5);
  size_t globalSize2D[] = {vDim / 4, vDim, 0};
  matAdd->configKernel(2, globalSize2D, nullptr, bufferA, bufferB, bufferC, &vDim, &vDim);

  for (int n = 0; n < 10; ++n) {
	float splitFactor = 1.0f / (float)pow2(n);
	splitKernelRun(vecAdd, exeQ, splitFactor);
  }

  memset(bufferC->hostPtr, 0, nByte);

  for (int n = 0; n < 10; ++n) {
	float splitFactor = 1.0f / (float)pow2(n);
	splitKernelRun(matAdd, exeQ, splitFactor);
  }

  delete edcl;
  return 0;
}

void splitKernelRun(EDKernel kernel, EDQueue exeQ, float splitFactor) {
  double timer;
  int numLoop = 1;
  edcl->confirmExeEnv(exeQ, kernel);
//  printGlobalSize(kernel);
  size_t gw[] = {static_cast<size_t>(ceil(kernel->globalSize[0] * splitFactor)),
				 kernel->globalSize[1],
				 kernel->globalSize[2]};
  kernel->setGlobalSize(gw);
  int numSplits = ceil(1 / splitFactor);
//  printGlobalSize(kernel);
  flushed_printf("Running Kernel %s, with splitFactor %f, ", kernel->name.c_str(), splitFactor);
  timer = getCurrentTime();
  for (int loopCnt = 0; loopCnt < numLoop; ++loopCnt) {
	for (int i = 0; i < numSplits; ++i) {
	  size_t k_offset[] = {i * gw[0], 0, 0};
	  kernel->setOffset(k_offset);
	  edcl->executeKernelBufferBlocked(exeQ, kernel, 0, nullptr, &kernel->kernelEvent);
	  kernel->waitForKernelEvent(exeQ);
	  //	flushed_printf("i = %d, Buffer C: \n", i);
	  //	printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
	}
  }
  timer = (getCurrentTime() - timer) / numLoop;
  flushed_printf("runtime: %lf\n", timer);
}

__unused void printGlobalSize(EDKernel kernel) {
  flushed_printf("Kernel %s Global Size: ", kernel->name.c_str());
  for (int i = 0; i < 3; ++i) flushed_printf("%d ", kernel->globalSize[i]);
  flushed_printf("\n");
}
