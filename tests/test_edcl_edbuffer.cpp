#pragma ide diagnostic ignored "DanglingPointers"
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include <unistd.h>
#include "EDCL.h"
#include "utils.h"

void EDBuffer_test(EDCL *edcl);

int main() {
  EDCL *edcl = new EDCL();
  EDBuffer_test(edcl);
  return 0;
}

void EDBuffer_test(EDCL *edcl) {
  size_t dataSize = 512 * 1024 * 1024;
  flushed_printf("EDBuffer cross platform accessing test, data size = %.2lfMB \n", (double)(dataSize >> 20));
  size_t numData = dataSize / sizeof(float);
  auto *reference = static_cast<float *>(malloc(dataSize));
  auto *cpuPulled = static_cast<float *>(malloc(dataSize));
  auto *gpuPulled = static_cast<float *>(malloc(dataSize));
  EDBuffer ed_buffer = edcl->createBuffer(dataSize);
  int numLoop = 5;
  for (int kI = 0; kI < numLoop; ++kI) {
	// copy, ref -> hostPtr
	// Init
	randArrayGenerator<float>(randNumberGeneratorCXX11(0, 100),
							  randNumberGeneratorCXX11(101, 1000),
							  reference,
							  numData);
	memset(cpuPulled, 0, dataSize);
	memset(gpuPulled, 0, dataSize);
	flushed_printf("ref -> hostPtr:");
	memcpy(ed_buffer->hostPtr, reference, dataSize);
	edcl->cpuRT_->pullBuffer(ed_buffer->CPUBuffer->mem, 0, dataSize, cpuPulled);
	edcl->gpuRT_->pullBuffer(ed_buffer->GPUzBuffer->mem, 0, dataSize, gpuPulled);
	flushed_printf(" memcmp(cpuPulled, reference):\t%d", memcmp(reference, cpuPulled, dataSize));
	flushed_printf(" memcmp(gpuPulled, reference):\t%d\n", memcmp(reference, gpuPulled, dataSize));

	// copy, ref -> CPUBuffer
	// Init
	randArrayGenerator<float>(randNumberGeneratorCXX11(0, 100),
							  randNumberGeneratorCXX11(101, 1000),
							  reference,
							  numData);
	memset(cpuPulled, 0, dataSize);
	memset(gpuPulled, 0, dataSize);
	flushed_printf("ref -> CPUBuffer:");
	edcl->cpuRT_->pushBuffer(ed_buffer->CPUBuffer, 0, dataSize, reference);
	edcl->cpuRT_->pullBuffer(ed_buffer->CPUBuffer->mem, 0, dataSize, cpuPulled);
	edcl->gpuRT_->pullBuffer(ed_buffer->GPUzBuffer->mem, 0, dataSize, gpuPulled);
	flushed_printf(" memcmp(hostPtr, reference):\t%d", memcmp(reference, ed_buffer->hostPtr, dataSize));
	flushed_printf(" memcmp(cpuPulled, reference):\t%d", memcmp(reference, cpuPulled, dataSize));
	flushed_printf(" memcmp(gpuPulled, reference):\t%d\n", memcmp(reference, gpuPulled, dataSize));

	// copy, ref -> GPUzBuffer
	// Init
	randArrayGenerator<float>(randNumberGeneratorCXX11(0, 100),
							  randNumberGeneratorCXX11(101, 1000),
							  reference,
							  numData);
	memset(cpuPulled, 0, dataSize);
	memset(gpuPulled, 0, dataSize);
	flushed_printf("ref -> GPUzBuffer:");
	edcl->gpuRT_->pushBuffer(ed_buffer->GPUzBuffer->mem, 0, dataSize, reference);
	edcl->cpuRT_->pullBuffer(ed_buffer->CPUBuffer, 0, dataSize, cpuPulled);
	edcl->gpuRT_->pullBuffer(ed_buffer->GPUzBuffer->mem, 0, dataSize, gpuPulled);
	flushed_printf(" memcmp(hostPtr, reference):\t%d", memcmp(reference, ed_buffer->hostPtr, dataSize));
	flushed_printf(" memcmp(cpuPulled, reference):\t%d", memcmp(reference, cpuPulled, dataSize));
	flushed_printf(" memcmp(gpuPulled, reference):\t%d\n", memcmp(reference, gpuPulled, dataSize));

  }

  free(cpuPulled);
  free(gpuPulled);
  free(reference);
  edcl->releaseEDBuffer(ed_buffer);
}
