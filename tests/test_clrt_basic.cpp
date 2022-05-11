#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include "testlib_cpu_vector_add.h"
#include "testlib_macros.h"
#include "CLRT.h"
#include "utils.h"
#include "CPU_utils.h"

const char *source =
#include "vecAdd.cl"

void clrt_vec_add_test(CLRT *clrt);
void clrt_mat_add_test(CLRT *clrt);

//#define PRINT_MAT

float *a;
float *b;
float *c;
float *c_cl;
#ifdef PRINT_MAT
unsigned int vDim = 4;
#else
unsigned int vDim = 8192;
#endif
unsigned int size = vDim * vDim;
size_t nByte = size * sizeof(float);
double timer;

int main() {
  flushed_printf("Initializing vector(%d)... ", size);
  a = (float *)malloc(nByte);
  b = (float *)malloc(nByte);
  TIME_CODE_ONCE(timer, randArrayGenerator<float>(0.001, 1, a, size); randArrayGenerator<float>(0.001, 1, b, size);)
//  for (int kI = 0; kI < size; ++kI) { a[kI] = ((float)kI + 1) / 100; }
//  memset(b, 0, nByte);
  c = (float *)malloc(nByte);
  c_cl = (float *)malloc(nByte);
  memset(c, 0, nByte);
  memset(c_cl, 0, nByte);
  flushed_printf("done. Time cost: %lf\n", timer);
  TIME_CODE_ONCE(timer, vecAddOptimized(a, b, c, size);)
  flushed_printf("Ground truth calculated, timer: %lf\n", timer);
#ifdef PRINT_MAT
  flushed_printf("A:\n");
  printMatrix(a, vDim, vDim);
  flushed_printf("B:\n");
  printMatrix(b, vDim, vDim);
  flushed_printf("C:\n");
  printMatrix(c, vDim, vDim);
#endif

  CLRT *cpuRT = CLRT::createCPURunTime();
  CLRT *gpuRT = CLRT::createGPURunTime();

  flushed_printf("CPU: %s\n", cpuRT->getOCLLibPath().c_str());
  std::vector<int> big = {4, 5, 6, 7};
  setThreadAffinity(big);
  printThreadAffinity(0);
  clrt_vec_add_test(cpuRT);
  clrt_mat_add_test(cpuRT);

  flushed_printf("GPU: %s\n", gpuRT->getOCLLibPath().c_str());
  std::vector<int> little = {1, 2, 3, 4};
  setThreadAffinity(little);
  printThreadAffinity(0);
  clrt_vec_add_test(gpuRT);
  clrt_mat_add_test(gpuRT);

  delete cpuRT;
  delete gpuRT;
  free(a);
  free(b);
  free(c);
  free(c_cl);
  return 0;
}

void clrt_vec_add_test(CLRT *clrt) {
  cl_event event;
  auto program = clrt->createProgramFromSource(&source);
  auto kernel_vacAdd = clrt->createKernel(program, "vecAdd", 4);
  cl_command_queue exeQ = clrt->createCLExecutionQueue(CL_QUEUE_PROFILING_ENABLE)->command_queue;
  auto bufferA = clrt->createBuffer(nByte, a);
  auto bufferB = clrt->createBuffer(nByte, nullptr);
  clrt->pushBuffer(bufferB, 0, nByte, b);
  auto bufferC = clrt->createBuffer(nByte, nullptr);
  auto arg_test_kernel = clrt->createKernel(program, "vecAdd", 4);
  CLRT_SET_ARGS(arg_test_kernel, bufferA, bufferB, bufferC, &size)
  clrt->setKernelArgWithLists(kernel_vacAdd,
							  arg_test_kernel->argTypeList,
							  arg_test_kernel->argPtrList,
							  arg_test_kernel->argSizeList,
							  arg_test_kernel->numArgs);
  size_t globalSize = size;
  printf("  OpenCL vecAdd:\n");
  // dry run
  clrt->execKernel(exeQ, kernel_vacAdd, 1,
				   nullptr, &globalSize, nullptr,
				   0, nullptr, &event);
  clWaitForEvents(1, &event, clrt->oclLibHandle_);
  clFinish(exeQ);
  clrt->releaseEvent(event);
  clrt->releaseEvent(event);
  // dry run end
  TIME_CODE_ONCE
  (timer,
   clrt->execKernel(exeQ, kernel_vacAdd, 1,
					nullptr, &globalSize, nullptr,
					0, nullptr, &event);
	   clWaitForEvents(1, &event, clrt->oclLibHandle_);
  )
  clFinish(exeQ);
  clrt->pullBuffer(bufferC, 0, size * sizeof(float), c_cl);
  flushed_printf("\tdiff = %f, memcmp = %d, ", calVecDiff(c, c_cl, size), memcmp(c, c_cl, nByte));
  flushed_printf("%lf sec, event time: %lf\n", timer, clrt->getExeTime(&event));

  flushed_printf("  Cleaning... ");
  flushed_printf("kernel_vecAdd ");
  clrt->releaseCLRT_Kernel(kernel_vacAdd);
  flushed_printf("vecA ");
  clrt->releaseCLRT_Buffer(bufferA);
  flushed_printf("bufferB ");
  clrt->releaseCLRT_Buffer(bufferB);
  flushed_printf("bufferC ");
  clrt->releaseCLRT_Buffer(bufferC);
  memset(c_cl, 0, nByte);
  flushed_printf("Done! \n\n");
}

void clrt_mat_add_test(CLRT *clrt) {
  cl_event event;
  CLRT_Program program = clrt->createProgramFromSource(&source);
  CLRT_Kernel kernel_matAdd = clrt->createKernel(program, "matAdd", 5);
  cl_command_queue exeQ = clrt->createCLExecutionProfilingQueue()->command_queue;
  CLRT_zBuffer bufferA = clrt->createZBuffer(nByte, a);
  CLRT_zBuffer bufferB = clrt->createZBuffer(nByte, nullptr);
  memcpy(bufferB->hostPtr, b, nByte);
  CLRT_zBuffer bufferC = clrt->createZBuffer(nByte, nullptr);

//  printMatrix(FLOAT_PTR(bufferA), vDim, vDim);
//  printMatrix(FLOAT_PTR(bufferB), vDim, vDim);
//  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);

  CLRT_SET_ARGS(kernel_matAdd, bufferA, bufferB, bufferC, &vDim, &vDim)
  size_t gw[] = {vDim / 4, vDim};
  printf("  OpenCL matAdd:\n");
  // dry run
  clrt->unmapZBuffer(bufferC);
  clrt->execKernel(exeQ, kernel_matAdd, 2,
				   nullptr, gw, nullptr,
				   0, nullptr, &event);
  clWaitForEvents(1, &event, clrt->oclLibHandle_);
  clFinish(exeQ);
  clrt->mapZBuffer(bufferC);
//  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
  memset(FLOAT_PTR(bufferC), 0, nByte);
  // dry run end
  clrt->unmapZBuffer(bufferC);
  TIME_CODE_ONCE
  (timer,
   clrt->execKernel(exeQ, kernel_matAdd, 2,
					nullptr, gw, nullptr,
					0, nullptr, &event);
	   clWaitForEvents(1, &event, clrt->oclLibHandle_);
  )
  clFinish(exeQ);
  clrt->mapZBuffer(bufferC);
//  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
#ifdef PRINT_MAT
  flushed_printf("A:\n");
  printMatrix((float *)vecA->hostPtr, vDim, vDim);
  flushed_printf("B:\n");
  printMatrix((float *)bufferB->hostPtr, vDim, vDim);
  flushed_printf("C:\n");
  printMatrix((float *)bufferC->hostPtr, vDim, vDim);
#endif
  flushed_printf("\tdiff = %f, memcmp = %d, ",
				 calVecDiff((float *)bufferC->hostPtr, c, size),
				 memcmp(bufferC->hostPtr, c, nByte));
  flushed_printf("%lf sec, event time: %lf\n", timer, clrt->getExeTime(&event));

  flushed_printf("  Cleaning... ");
  flushed_printf("kernel_matAdd ");
  clrt->releaseCLRT_Kernel(kernel_matAdd);
  flushed_printf("vecA ");
  clrt->releaseCLRT_zBuffer(bufferA);
  flushed_printf("bufferB ");
  clrt->releaseCLRT_zBuffer(bufferB);
  flushed_printf("bufferC ");
  clrt->releaseCLRT_zBuffer(bufferC);
  memset(c_cl, 0, nByte);
  flushed_printf("Done! \n\n");
}
