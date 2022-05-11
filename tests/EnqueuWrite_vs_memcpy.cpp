#include <cstdlib>
#include "CLRT.h"
#include "utils.h"
#include <cmath>

#define MB 1024*1024

#define TIME_CODE_LOOP(TIMER, LOOP_NUM, CODE_BLOCK) \
{ \
    TIMER=getCurrentTime(); \
    for (int loopCnt = 0; loopCnt < LOOP_NUM; ++loopCnt) \
        {CODE_BLOCK} \
    TIMER=(getCurrentTime()-TIMER)/LOOP_NUM; \
}

void creationTest(CLRT *clrt)__attribute__((optnone)) {
  flushed_printf("  Creation & free test(clCreateBuffer vs calloc )... \n");
  unsigned int loopNum;
  for (int i = 1; i < 10; ++i) {
	loopNum = i > 8 ? 100 : 1000;
	size_t dataSize = pow2(i) * MB;
	double cl_timer;
	TIME_CODE_LOOP(cl_timer, loopNum,
				   CLRT_zBuffer zBuffer = clrt->createZBuffer(dataSize, nullptr);
					   clrt->releaseCLRT_zBuffer(zBuffer);
	)
	double c_timer;
	TIME_CODE_LOOP(c_timer, loopNum,
				   char *ptr = (char *)calloc(dataSize, 1);
					   free(ptr);)
	flushed_printf("\t%d MB,\t%lf, %lf, %lf, %s\n",
				   pow2(i),
				   cl_timer,
				   c_timer,
				   cl_timer / c_timer,
				   cl_timer > c_timer ? "calloc" : "clCreateBuffer");
  }
}

void pushTest(CLRT *clrt, void *data_ptr) {
  flushed_printf("  Push test(clEnqueueWrite vs memcpy)... \n");
  unsigned int loopNum = 100;
  for (int i = 1; i < 10; ++i) {
	size_t dataSize = pow2(i) * MB;
	CLRT_zBuffer zBuffer = clrt->createZBuffer(dataSize, nullptr);
	double cl_timer = 0;
	TIME_CODE_LOOP(cl_timer, loopNum,
				   clrt->pushBuffer(zBuffer->mem, 0, dataSize, data_ptr);)
	double c_timer = 0;
	TIME_CODE_LOOP(c_timer, loopNum,
				   memcpy(zBuffer->hostPtr, data_ptr, dataSize);)
	flushed_printf("\t%d MB,\t%lf, %lf, %lf, %s\n",
				   pow2(i),
				   cl_timer,
				   c_timer,
				   cl_timer / c_timer,
				   cl_timer > c_timer ? "memcpy_w" : "clEnqueueWrite");
	clrt->releaseCLRT_zBuffer(zBuffer);
  }
}

void pullTest(CLRT *clrt, void *data_ptr) {
  flushed_printf("  Pull test(clEnqueueRead vs memcpy)... \n");
  unsigned int loopNum = 100;
  for (int i = 1; i < 10; ++i) {
	size_t dataSize = pow2(i) * MB;
	CLRT_zBuffer zBuffer = clrt->createZBuffer(dataSize, nullptr);
	double cl_timer = 0;
	TIME_CODE_LOOP(cl_timer, loopNum,
				   clrt->pullBuffer(zBuffer->mem, 0, dataSize, data_ptr);)
	double c_timer = 0;
	TIME_CODE_LOOP(c_timer, loopNum,
				   memcpy(data_ptr, zBuffer->hostPtr, dataSize);)
	flushed_printf("\t%d MB,\t%lf, %lf, %lf, %s\n",
				   pow2(i),
				   cl_timer,
				   c_timer,
				   cl_timer / c_timer,
				   cl_timer > c_timer ? "memcpy_r" : "clEnqueueWrite");
	clrt->releaseCLRT_zBuffer(zBuffer);
  }
}

void accessTest(CLRT *clrt) {
  flushed_printf("  Access test(clCreateBuffer(CL_MEM_ALLOC_HOST_PTR) host_ptr vs calloc ptr)... \n");
  unsigned int loopNum = 10;
  for (int i = 1; i < 10; ++i) {
	size_t dataSize = pow2(i) * MB;
	CLRT_zBuffer zBuffer = clrt->createZBuffer(dataSize, nullptr);
	char *host_ptr = (char *)zBuffer->hostPtr;
	char *ptr = (char *)calloc(dataSize, sizeof(char));
	double cl_timer = 0;
	TIME_CODE_LOOP(cl_timer, loopNum,
				   for (int j = 0; j < dataSize; ++j) host_ptr[j]++;)
	double c_timer = 0;
	TIME_CODE_LOOP(c_timer, loopNum,
				   for (int j = 0; j < dataSize; ++j) ptr[j]++;)
	flushed_printf("\t%d MB,\t%lf, %lf, %lf, %s\n",
				   pow2(i),
				   cl_timer,
				   c_timer,
				   cl_timer / c_timer,
				   cl_timer > c_timer ? "memcpy_r" : "clEnqueueWrite");
	clrt->releaseCLRT_zBuffer(zBuffer);
  }
}

int main() {
  void *host_ptr = malloc(512 * MB);
  CLRT *CPUrt = CLRT::createCPURunTime();
  flushed_printf("POCL: \n");
  creationTest(CPUrt);
  pushTest(CPUrt, host_ptr);
  pullTest(CPUrt, host_ptr);
  accessTest(CPUrt);
  CLRT *GPUrt = CLRT::createGPURunTime();
  flushed_printf("GPU: \n");
  creationTest(GPUrt);
  pushTest(GPUrt, host_ptr);
  pullTest(GPUrt, host_ptr);
  accessTest(GPUrt);
  delete CPUrt;
  delete GPUrt;
}
