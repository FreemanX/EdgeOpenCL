#include <tar.h>
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
// cross platform buffer prototype
#include "utils.h"
#include "CLRT.h"

//  C_BUF: cross platform buffer. GPU Z_BUF + CPU BUF
#define C_BUF 0x7675625746C75
typedef struct CBufferTypeA_ {
  __unused long SIGNATURE = C_BUF;
  CLRT_Buffer CPU_Buffer{};
  CLRT_zBuffer GPU_zBuffer{};
  void *host_ptr{};
  size_t sizeByte{};
} CBufferTypeA_;
typedef struct CBufferTypeA_ *CLRT_cBufferA;

typedef struct CBufferTypeB_ {
  __unused long SIGNATURE = C_BUF;
  CLRT_zBuffer CPU_zBuffer{};
  CLRT_Buffer GPU_Buffer{};
  void *host_ptr{};
  size_t sizeByte{};
} CBufferTypeB_;

typedef struct CBufferTypeA_ *CLRT_cBufferA;
typedef struct CBufferTypeB_ *CLRT_cBufferB;

CLRT_cBufferA createCBufferA(CLRT *CPUrt, CLRT *GPUrt, size_t nByte) {
  auto c_buffer = (CLRT_cBufferA)calloc(1, sizeof(CBufferTypeA_));
  c_buffer->SIGNATURE = C_BUF;
  c_buffer->GPU_zBuffer = GPUrt->createZBuffer(nByte, nullptr);
  c_buffer->host_ptr = c_buffer->GPU_zBuffer->hostPtr;
  c_buffer->CPU_Buffer = CPUrt->createBuffer(CL_MEM_USE_HOST_PTR, nByte, c_buffer->host_ptr);
  c_buffer->sizeByte = nByte;
  return c_buffer;
}

void releaseCBufferA(CLRT *CPUrt, CLRT *GPUrt, CLRT_cBufferA c_buffer) {
  if (c_buffer->sizeByte) {
	CPUrt->releaseCLRT_Buffer(c_buffer->CPU_Buffer);
	GPUrt->releaseCLRT_zBuffer(c_buffer->GPU_zBuffer);
	c_buffer->sizeByte = 0;
	free(c_buffer);
  }
}

CLRT_cBufferB createCBufferB(CLRT *CPUrt, CLRT *GPUrt, size_t nByte) {
  auto c_buffer = (CLRT_cBufferB)calloc(1, sizeof(CBufferTypeB_));
  c_buffer->SIGNATURE = C_BUF;
  c_buffer->CPU_zBuffer = CPUrt->createZBuffer(nByte, nullptr);
  c_buffer->host_ptr = c_buffer->CPU_zBuffer->hostPtr;
  c_buffer->GPU_Buffer = GPUrt->createBuffer(CL_MEM_USE_HOST_PTR, nByte, c_buffer->host_ptr);
  c_buffer->sizeByte = nByte;
  return c_buffer;
}

void releaseCBufferB(CLRT *CPUrt, CLRT *GPUrt, CLRT_cBufferB c_buffer) {
  if (c_buffer->sizeByte) {
	GPUrt->releaseCLRT_Buffer(c_buffer->GPU_Buffer);
	CPUrt->releaseCLRT_zBuffer(c_buffer->CPU_zBuffer);
	c_buffer->sizeByte = 0;
	free(c_buffer);
  }
}

int main() {
  CLRT *CPUrt = CLRT::createCPURunTime();
  CLRT *GPUrt = CLRT::createGPURunTime();
  size_t dataSize = 512 * 1024 * 1024;
  size_t numData = dataSize / sizeof(float);
  auto *reference = (float *)malloc(dataSize);
  randArrayGenerator<float>(0.001, 1, reference, numData);
//  for (int i = 0; i < numData; ++i) reference[i] = (float)i;

  auto *poclPulled = (float *)calloc(numData, sizeof(float));
  auto *GPUPulled = (float *)calloc(numData, sizeof(float));

  // Operate on c_buffer's host_ptr should changed both cl_mem for CPU and GPU
  flushed_printf("TypeA, CPU Buffer, GPU zBuffer\n");
  CLRT_cBufferA c_bufferA = createCBufferA(CPUrt, GPUrt, dataSize);
  flushed_printf("Compare poclPulled before pull: %d\n", memcmp(poclPulled, reference, dataSize));
  flushed_printf("Compare GPUPulled  before pull: %d\n", memcmp(GPUPulled, reference, dataSize));
  memcpy(c_bufferA->host_ptr, reference, dataSize);
  CPUrt->pullBuffer(c_bufferA->CPU_Buffer->mem, 0, dataSize, poclPulled);
  flushed_printf("Compare poclPulled: %d\n", memcmp(poclPulled, reference, dataSize));
  GPUrt->pullBuffer(c_bufferA->GPU_zBuffer->mem, 0, dataSize, GPUPulled);
  flushed_printf("Compare GPUPulled: %d\n", memcmp(GPUPulled, reference, dataSize));
  releaseCBufferA(CPUrt, GPUrt, c_bufferA);

  memset(poclPulled, 0, dataSize);
  memset(GPUPulled, 0, dataSize);
  flushed_printf("TypeB, CPU zBuffer, GPU Buffer\n");
  CLRT_cBufferB c_bufferB = createCBufferB(CPUrt, GPUrt, dataSize);
  flushed_printf("Compare poclPulled before pull: %d\n", memcmp(poclPulled, reference, dataSize));
  flushed_printf("Compare GPUPulled  before pull: %d\n", memcmp(GPUPulled, reference, dataSize));
  memcpy(c_bufferB->host_ptr, reference, dataSize);
  CPUrt->pullBuffer(c_bufferB->CPU_zBuffer->mem, 0, dataSize, poclPulled);
  flushed_printf("Compare poclPulled: %d\n", memcmp(poclPulled, reference, dataSize));
  GPUrt->pullBuffer(c_bufferB->GPU_Buffer->mem, 0, dataSize, GPUPulled);
  flushed_printf("Compare GPUPulled: %d\n", memcmp(GPUPulled, reference, dataSize));
  releaseCBufferB(CPUrt, GPUrt, c_bufferB);
  flushed_printf("Done !\n");
  return 0;
}
