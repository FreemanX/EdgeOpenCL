#include <unistd.h>
#include "EDCL.h"

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
EDQueue exeQ;
EDProgram program;
EDKernel vecAdd;
size_t globalVecAdd[1] = {size};
size_t localVecAdd[1] = {1024};
EDKernel matAdd;
size_t globalMatAdd[2] = {vDim / 4, vDim};
size_t localMatAdd[2] = {32, 32};

void init() {
  edcl = new EDCL();
  bufferA = edcl->createBuffer(nByte);
  bufferB = edcl->createBuffer(nByte);
//  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferA), size);
//  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferB), size);
//  for (int i = 0; i < size; ++i) { FLOAT_PTR(bufferA)[i] = (float)(i + 1) / 100; }
  for (int i = 0; i < size; ++i) { FLOAT_PTR(bufferA)[i] = 1; }
  memset(bufferB->hostPtr, 0, nByte);
  bufferC = edcl->createBuffer(nByte);
  memset(bufferC->hostPtr, 0, nByte);
//  exeQ = edcl->createDeviceCmdQueueProfilingEnabled(CPUQueue);
//  exeQ = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  exeQ = edcl->createSubDeviceExeQueue(2, 1);
  program = edcl->createProgram(&source);
  vecAdd = edcl->createKernel(program, "vecAdd", 4);
  vecAdd->configKernel(1, globalVecAdd, localVecAdd, bufferA, bufferB, bufferC, &size);
  matAdd = edcl->createKernel(program, "matAdd", 5);
  matAdd->configKernel(2, globalMatAdd, localMatAdd, bufferA, bufferB, bufferC, &vDim, &vDim);
}

//AtomicKernelSet slicingKernel2D(EDKernel k, uint sliceFactor) {
//  printf("\033[1;34m");
//  if (k->localSize == nullptr) return nullptr;
//  uint totalNumSlices = pow2(sliceFactor);
//  size_t numWGs[2]{k->globalSize[0] / k->localSize[0], k->globalSize[1] / k->localSize[1]};
//  if (numWGs[0] * numWGs[1] < totalNumSlices) {
//	flushed_printf("Slice factor %d exceeds the kernel global size. ", sliceFactor);
//	flushed_printf("Number of WGs({%d, %d}):%d. Required number of slices: %d.\n",
//				   numWGs[0], numWGs[1], numWGs[0] * numWGs[1], totalNumSlices);
//	return nullptr;
//  }
//  flushed_printf("Number of WGs({%d, %d}):%d. Required number of slices: %d.\n",
//				 numWGs[0],
//				 numWGs[1],
//				 numWGs[0] * numWGs[1],
//				 totalNumSlices);
//  uint numSlices[2]; // number of slices on each dim
//  for (int i = 0; i < sliceFactor || sliceFactor == 0; i++) {
//	numSlices[0] = pow2(sliceFactor - i);
//	numSlices[1] = pow2(i);
//	if (numWGs[0] >= numSlices[0] && numWGs[1] >= numSlices[1]) { break; }
//  }
//  flushed_printf("numSlices[%d, %d]\n", numSlices[0], numSlices[1]);
//  auto *sliceGlobalSize = new size_t[2]; // each slice will have this many WGs in dim0 and dim1
//  sliceGlobalSize[0] = (numWGs[0] / numSlices[0]) * k->localSize[0]; // number of WGs per slice * WG size
//  sliceGlobalSize[1] = (numWGs[1] / numSlices[1]) * k->localSize[1];
//  AtomicKernelSet aks = edcl->createAtomicKernelSet(totalNumSlices);
//  for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//	for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//	  EDKernel kernelSlice = edcl->createKernelCopy(k);
//	  kernelSlice->globalSize = sliceGlobalSize;
//	  kernelSlice->offset = new size_t[2]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1};
//	  aks->addKernel(kernelSlice);
//	}
//  }
//  return aks;
//}
//
//AtomicKernelSet slicingKernelPow2(EDKernel k, uint sliceFactor) {
//  printf("\033[1;34m%s: ", __func__);
//  if (k->localSize == nullptr) {
//	PRINT_RED
//	flushed_printf("Kernel %s local work size not set!\n", k->name.data());
//	PRINT_DEFAULT
//	return nullptr;
//  }
//  uint kDim = k->workDim; // kernel NDRange dimensions
//  for (int kI = 0; kI < kDim; ++kI) {
//	if (k->globalSize[kI] % k->localSize[kI] != 0) {
//	  PRINT_RED
//	  flushed_printf("Kernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n",
//					 k->name.data(),
//					 kI, k->globalSize[kI], k->localSize[kI]);
//	  PRINT_DEFAULT
//	  return nullptr;
//	}
//  }
//  uint requiredNumSlices = pow2(sliceFactor);
//  auto *numWGs = new size_t[kDim];
//  uint maxNumSlices = 1;
//  for (int kI = 0; kI < kDim; ++kI) {
//	numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
//	maxNumSlices *= numWGs[kI];
//  }
//  if (maxNumSlices < requiredNumSlices) {
//	PRINT_RED
//	flushed_printf("Slice factor %d exceeds the kernel global size. ", sliceFactor);
//	flushed_printf("Number of WGs({");
//	for (int kI = 0; kI < kDim; ++kI)
//	  flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
//	flushed_printf("}):%d. Required number of slices: %d.\n", maxNumSlices, requiredNumSlices);
//	PRINT_DEFAULT
//	return nullptr;
//  }
////  flushed_printf("Max Number of WGs({");
////  for (int kI = 0; kI < kDim; ++kI)
////	flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
////  flushed_printf("}):%d. Required number of slices: %d.\n", maxNumSlices, requiredNumSlices);
//
//  auto *numSlices = new size_t[kDim];
//  for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
//  int remainSliceFactor = sliceFactor;
//  for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
//	for (int i = 0; i < remainSliceFactor; ++i) {
//	  uint dimSliceFactor = remainSliceFactor - i;
//	  if (numWGs[kI] >= pow2(dimSliceFactor)) {
//		numSlices[kI] = pow2(dimSliceFactor);
//		remainSliceFactor -= dimSliceFactor;
//		break;
//	  }
//	}
////	flushed_printf("remainSliceFactor: %d\n", remainSliceFactor);
//  }
////  flushed_printf("numSlices[");
////  for (int kI = 0; kI < kDim; ++kI)
////	flushed_printf("%d%s", numSlices[kI], kI + 1 == kDim ? "" : ", ");
////  flushed_printf("]\n", maxNumSlices, requiredNumSlices);
//
//  auto *sliceGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
////  flushed_printf("sliceGlobalSize[ ");
//  for (int kI = 0; kI < kDim; ++kI) {
//	sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
////	flushed_printf("%d ", sliceGlobalSize[kI]);
//  }
////  flushed_printf("]\n");
//
//  AtomicKernelSet aks = edcl->createAtomicKernelSet(requiredNumSlices);
//  if (kDim == 1) {
//	for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//	  EDKernel kernelSlice = edcl->createKernelCopy(k);
//	  kernelSlice->globalSize = sliceGlobalSize;
//	  kernelSlice->offset = new size_t[1]{sliceGlobalSize[0] * dim0};
//	  if (sliceFactor != 0 && dim0 == numSlices[0] - 1) {
//		auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//		remainGlobalSize[0] = sliceGlobalSize[0];
//		if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//		  remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
//		}
//		kernelSlice->globalSize = remainGlobalSize;
//	  }
//	  aks->addKernel(kernelSlice);
//	}
//  } else if (kDim == 2) {
//	for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//	  for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		EDKernel kernelSlice = edcl->createKernelCopy(k);
//		kernelSlice->globalSize = sliceGlobalSize;
//		kernelSlice->offset = new size_t[kDim]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1};
//		if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1)) { // if it's the edge slice
//		  auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//		  remainGlobalSize[0] = sliceGlobalSize[0];
//		  remainGlobalSize[1] = sliceGlobalSize[1];
//		  if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
//		  }
//		  if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
//		  }
//		  kernelSlice->globalSize = remainGlobalSize;
////		  flushed_printf("Kernel %s globalSize[%d,%d], kernelSlice->offset[%d,%d], ",
////						 k->name.data(), k->globalSize[0], k->globalSize[1],
////						 kernelSlice->offset[0], kernelSlice->offset[1]);
////		  flushed_printf("Slice[%d,%d] kernelSlice->globalSize[%d,%d]\n",
////						 dim0,
////						 dim1,
////						 kernelSlice->globalSize[0],
////						 kernelSlice->globalSize[0]);
//		}
//		aks->addKernel(kernelSlice);
//	  }
//	}
//  } else {
//	for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
//	  for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		  EDKernel kernelSlice = edcl->createKernelCopy(k);
//		  kernelSlice->globalSize = sliceGlobalSize;
//		  kernelSlice->offset =
//			  new size_t[3]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1, sliceGlobalSize[2] * dim2};
//		  if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
//			  || dim2 == numSlices[2] - 1)) { // if it's the edge slice
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = sliceGlobalSize[0];
//			remainGlobalSize[1] = sliceGlobalSize[1];
//			remainGlobalSize[2] = sliceGlobalSize[2];
//			if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			  remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
//			}
//			if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			  remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
//			}
//			if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
//			  remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
//			}
//			kernelSlice->globalSize = remainGlobalSize;
//		  }
//		  aks->addKernel(kernelSlice);
//		}
//	  }
//	}
//  }
//  delete[] numWGs;
//  delete[] numSlices;
//  printf("\033[0m");
//  return aks;
//}
//
//// create a extra slices if there's remaining work
//AtomicKernelSet slicingKernelPow2_extra(EDKernel k, uint sliceFactor) {
//  printf("\033[1;34m%s: ", __func__);
//  if (k->localSize == nullptr) {
//	PRINT_RED
//	flushed_printf("Kernel %s local work size not set!\n", k->name.data());
//	PRINT_DEFAULT
//	return nullptr;
//  }
//  uint kDim = k->workDim; // kernel NDRange dimensions
//  for (int kI = 0; kI < kDim; ++kI) {
//	if (k->globalSize[kI] % k->localSize[kI] != 0) {
//	  PRINT_RED
//	  flushed_printf("Kernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n",
//					 k->name.data(),
//					 kI, k->globalSize[kI], k->localSize[kI]);
//	  PRINT_DEFAULT
//	  return nullptr;
//	}
//  }
//  uint requiredNumSlices = pow2(sliceFactor);
//  auto *numWGs = new size_t[kDim];
//  uint maxNumSlices = 1;
//  for (int kI = 0; kI < kDim; ++kI) {
//	numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
//	maxNumSlices *= numWGs[kI];
//  }
//  if (maxNumSlices < requiredNumSlices) {
//	PRINT_RED
//	flushed_printf("Slice factor %d exceeds the kernel global size. ", sliceFactor);
//	flushed_printf("Number of WGs({");
//	for (int kI = 0; kI < kDim; ++kI)
//	  flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
//	flushed_printf("}):%d. Required number of slices: %d.\n", maxNumSlices, requiredNumSlices);
//	PRINT_DEFAULT
//	return nullptr;
//  }
//
//  auto *numSlices = new size_t[kDim];
//  for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
//  int remainSliceFactor = sliceFactor;
//  for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
//	for (int i = 0; i < remainSliceFactor; ++i) {
//	  uint dimSliceFactor = remainSliceFactor - i;
//	  if (numWGs[kI] >= pow2(dimSliceFactor)) {
//		numSlices[kI] = pow2(dimSliceFactor);
//		remainSliceFactor -= dimSliceFactor;
//		break;
//	  }
//	}
//  }
////  flushed_printf("numSlices[");
////  for (int kI = 0; kI < kDim; ++kI)
////	flushed_printf("%d%s", numSlices[kI], kI + 1 == kDim ? "" : ", ");
////  flushed_printf("]\n", maxNumSlices, requiredNumSlices);
//
//  auto *sliceGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//  for (int kI = 0; kI < kDim; ++kI) {
//	sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
//  }
//
//  AtomicKernelSet aks = edcl->createAtomicKernelSet(requiredNumSlices);
//  if (kDim == 1) {
//	for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//	  EDKernel kernelSlice = edcl->createKernelCopy(k);
//	  kernelSlice->globalSize = sliceGlobalSize;
//	  kernelSlice->offset = new size_t[1]{sliceGlobalSize[0] * dim0};
//	  aks->addKernel(kernelSlice);
//	  if (sliceFactor != 0) {
//		if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//		  EDKernel extraSlice = edcl->createKernelCopy(k);
//		  extraSlice->offset = new size_t[1]{kernelSlice->offset[0] + kernelSlice->globalSize[0]};
//		  extraSlice->globalSize = new size_t[1]{k->globalSize[0] - extraSlice->offset[0]};
//		  aks->addKernel(extraSlice);
//		}
//	  }
//	}
//  } else if (kDim == 2) {
//	for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//	  for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		EDKernel kernelSlice = edcl->createKernelCopy(k);
//		kernelSlice->globalSize = sliceGlobalSize;
//		kernelSlice->offset = new size_t[kDim]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1};
//		aks->addKernel(kernelSlice);
////		/*
//		if (sliceFactor != 0) { // if it's the edge slice
//		  if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
//			remainGlobalSize[1] = sliceGlobalSize[1];
//			auto *remainOffset = (size_t *)calloc(kDim, sizeof(size_t));
//			remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
//			remainOffset[1] = kernelSlice->offset[1];
//
//			EDKernel extraSlice = edcl->createKernelCopy(k);
//			extraSlice->offset = remainOffset;
//			extraSlice->globalSize = remainGlobalSize;
//			aks->addKernel(extraSlice);
//		  }
//		  if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = sliceGlobalSize[0];
//			remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
//			auto *remainOffset = (size_t *)calloc(kDim, sizeof(size_t));
//			remainOffset[0] = kernelSlice->offset[0];
//			remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
//
//			EDKernel extraSlice = edcl->createKernelCopy(k);
//			extraSlice->offset = remainOffset;
//			extraSlice->globalSize = remainGlobalSize;
//			aks->addKernel(extraSlice);
//		  }
//		  if ((dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1])
//			  && dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
//			remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
//			auto *remainOffset = (size_t *)calloc(kDim, sizeof(size_t));
//			remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
//			remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
//
//			EDKernel extraSlice = edcl->createKernelCopy(k);
//			extraSlice->offset = remainOffset;
//			extraSlice->globalSize = remainGlobalSize;
//			aks->addKernel(extraSlice);
//		  }
//		}
////		 */
//	  }
//	}
//  } else {
//	for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
//	  for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		  EDKernel kernelSlice = edcl->createKernelCopy(k);
//		  kernelSlice->globalSize = sliceGlobalSize;
//		  kernelSlice->offset =
//			  new size_t[3]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1, sliceGlobalSize[2] * dim2};
//		  if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
//			  || dim2 == numSlices[2] - 1)) { // if it's the edge slice
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = sliceGlobalSize[0];
//			remainGlobalSize[1] = sliceGlobalSize[1];
//			remainGlobalSize[2] = sliceGlobalSize[2];
//			if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			  remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
//			}
//			if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			  remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
//			}
//			if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
//			  remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
//			}
//			kernelSlice->globalSize = remainGlobalSize;
//		  }
//		  aks->addKernel(kernelSlice);
//		}
//	  }
//	}
//  }
//  delete[] numWGs;
//  delete[] numSlices;
//  printf("\033[0m");
//  return aks;
//}
//
//// create equal size extra slices(try) if there's remaining work
//AtomicKernelSet slicingKernelPow2_extra_equal(EDKernel k, uint sliceFactor) {
//  printf("\033[1;34m%s: ", __func__);
//  if (k->localSize == nullptr) {
//	flushed_printf("\033[1;31mKernel %s local work size not set!\n\033[0m", k->name.data());
//	return nullptr;
//  }
//  uint kDim = k->workDim; // kernel NDRange dimensions
//  for (int kI = 0; kI < kDim; ++kI) {
//	if (k->globalSize[kI] % k->localSize[kI] != 0) {
//	  flushed_printf(
//		  "\033[1;31mKernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n\033[0m",
//		  k->name.data(),
//		  kI,
//		  k->globalSize[kI],
//		  k->localSize[kI]);
//	  return nullptr;
//	}
//  }
//  uint requiredNumSlices = pow2(sliceFactor);
//  auto *numWGs = new size_t[kDim];
//  uint maxNumSlices = 1;
//  for (int kI = 0; kI < kDim; ++kI) {
//	numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
//	maxNumSlices *= numWGs[kI];
//  }
//  if (maxNumSlices < requiredNumSlices) {
//	flushed_printf("\033[1;31mSlice factor %d exceeds the kernel global size. ", sliceFactor);
//	flushed_printf("Number of WGs({");
//	for (int kI = 0; kI < kDim; ++kI)
//	  flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
//	flushed_printf("}):%d. Required number of slices: %d.\n\033[0m", maxNumSlices, requiredNumSlices);
//	return nullptr;
//  }
//
//  auto *numSlices = new size_t[kDim];
//  for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
//  int remainSliceFactor = sliceFactor;
//  for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
//	for (int i = 0; i < remainSliceFactor; ++i) {
//	  uint dimSliceFactor = remainSliceFactor - i;
//	  if (numWGs[kI] >= pow2(dimSliceFactor)) {
//		numSlices[kI] = pow2(dimSliceFactor);
//		remainSliceFactor -= dimSliceFactor;
//		break;
//	  }
//	}
//  }
////  flushed_printf("numSlices[");
////  for (int kI = 0; kI < kDim; ++kI)
////	flushed_printf("%d%s", numSlices[kI], kI + 1 == kDim ? "" : ", ");
////  flushed_printf("]\n", maxNumSlices, requiredNumSlices);
//
//  auto *sliceGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//  for (int kI = 0; kI < kDim; ++kI) {
//	sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
//  }
//
//  AtomicKernelSet aks = edcl->createAtomicKernelSet(requiredNumSlices);
//  if (kDim == 1) {
//	for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//	  EDKernel kernelSlice = edcl->createKernelCopy(k);
//	  kernelSlice->globalSize = sliceGlobalSize;
//	  kernelSlice->offset = new size_t[1]{sliceGlobalSize[0] * dim0};
//	  aks->addKernel(kernelSlice);
//	  if (sliceFactor != 0) {
//		size_t remainWI0 = k->globalSize[0] - (kernelSlice->offset[0] + sliceGlobalSize[0]);
//		if (dim0 == numSlices[0] - 1 && remainWI0 > 0) {
//		  size_t remainChunk = ceil((float)remainWI0 / (float)sliceGlobalSize[0]);
//		  for (int i = 0; i < remainChunk; ++i) {
//			EDKernel extraSlice = edcl->createKernelCopy(k);
//			extraSlice->offset = new size_t[1]{kernelSlice->offset[0] + (i + 1) * sliceGlobalSize[0]};
//			extraSlice->globalSize = new size_t[1];
//			extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
//			aks->addKernel(extraSlice);
//		  }
//		}
//	  }
//	}
//  } else if (kDim == 2) {
//	for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//	  for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		EDKernel kernelSlice = edcl->createKernelCopy(k);
//		kernelSlice->globalSize = sliceGlobalSize;
//		kernelSlice->offset = new size_t[kDim]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1};
//		aks->addKernel(kernelSlice);
////		/*
//		if (sliceFactor != 0) { // if it's the edge slice
//		  if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			size_t remainGlobalSize[2];
//			remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
//			remainGlobalSize[1] = sliceGlobalSize[1];
//			size_t remainOffset[2];
//			remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
//			remainOffset[1] = kernelSlice->offset[1];
//
//			size_t remainChunk = ceil((float)remainGlobalSize[0] / (float)sliceGlobalSize[0]);
//			for (int i = 0; i < remainChunk; ++i) {
//			  EDKernel extraSlice = edcl->createKernelCopy(k);
//			  extraSlice->offset = new size_t[2]{remainOffset[0] + sliceGlobalSize[0] * i, remainOffset[1]};
//			  extraSlice->globalSize = new size_t[2];
//			  extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
//			  extraSlice->globalSize[1] = remainGlobalSize[1];
//			  aks->addKernel(extraSlice);
//			}
//
//		  }
//		  if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			size_t remainGlobalSize[2];
//			remainGlobalSize[0] = sliceGlobalSize[0];
//			remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
//			size_t remainOffset[2];
//			remainOffset[0] = kernelSlice->offset[0];
//			remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
//
//			size_t remainChunk = ceil((float)remainGlobalSize[1] / (float)sliceGlobalSize[1]);
//			for (int i = 0; i < remainChunk; ++i) {
//			  EDKernel extraSlice = edcl->createKernelCopy(k);
//			  extraSlice->offset = new size_t[2]{remainOffset[0], remainOffset[1] + sliceGlobalSize[1] * i};
//			  extraSlice->globalSize = new size_t[2];
//			  extraSlice->globalSize[0] = remainGlobalSize[0];
//			  extraSlice->globalSize[1] = MIN(sliceGlobalSize[1], k->globalSize[1] - extraSlice->offset[1]);
//			  aks->addKernel(extraSlice);
//			}
//		  }
//		  if ((dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1])
//			  && dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			size_t remainGlobalSize[2];
//			remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
//			remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
//			size_t remainOffset[2];
//			remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
//			remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
//
//			size_t remainChunk[2];
//			remainChunk[0] = ceil((float)remainGlobalSize[0] / (float)sliceGlobalSize[0]);
//			remainChunk[1] = ceil((float)remainGlobalSize[1] / (float)sliceGlobalSize[1]);
//			for (int i = 0; i < remainChunk[0]; ++i) {
//			  for (int j = 0; j < remainChunk[1]; ++j) {
//				EDKernel extraSlice = edcl->createKernelCopy(k);
//				extraSlice->offset =
//					new size_t[2]{remainOffset[0] + sliceGlobalSize[0] * i, remainOffset[1] + sliceGlobalSize[1] * j};
//				extraSlice->globalSize = new size_t[2];
//				extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
//				extraSlice->globalSize[1] = MIN(sliceGlobalSize[1], k->globalSize[1] - extraSlice->offset[1]);
//				aks->addKernel(extraSlice);
//			  }
//			}
//		  }
//		}
////		 */
//	  }
//	}
//  } else {
//	for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
//	  for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
//		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
//		  EDKernel kernelSlice = edcl->createKernelCopy(k);
//		  kernelSlice->globalSize = sliceGlobalSize;
//		  kernelSlice->offset =
//			  new size_t[3]{sliceGlobalSize[0] * dim0, sliceGlobalSize[1] * dim1, sliceGlobalSize[2] * dim2};
//		  if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
//			  || dim2 == numSlices[2] - 1)) { // if it's the edge slice
//			auto *remainGlobalSize = (size_t *)calloc(kDim, sizeof(size_t));
//			remainGlobalSize[0] = sliceGlobalSize[0];
//			remainGlobalSize[1] = sliceGlobalSize[1];
//			remainGlobalSize[2] = sliceGlobalSize[2];
//			if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
//			  remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
//			}
//			if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
//			  remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
//			}
//			if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
//			  remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
//			}
//			kernelSlice->globalSize = remainGlobalSize;
//		  }
//		  aks->addKernel(kernelSlice);
//		}
//	  }
//	}
//  }
//  delete[] numWGs;
//  delete[] numSlices;
//  printf("\033[0m");
//  return aks;
//}

void testSlicingDebug(EDKernel kernel, AtomicKernelSet sliceFunc(EDKernel, uint)) {
  uint sliceFactor = 0;
  auto aks = sliceFunc(kernel, sliceFactor);
  while (aks != nullptr) {
	for (int i = 0; i < aks->getNumKernels(); i++) {
	  EDKernel k = aks->kernels[i];
//	  flushed_printf("\033[2J"); // clear screen
	  flushed_printf("offset[ ");
	  for (int kI = 0; kI < k->workDim; ++kI) {
		flushed_printf("%d ", k->offset[kI]);
	  }
	  flushed_printf("], ");
	  flushed_printf("GlobalSize[ ");
	  for (int kI = 0; kI < k->workDim; ++kI) {
		flushed_printf("%d ", k->globalSize[kI]);
	  }
	  flushed_printf("]\n");
	  edcl->confirmExeEnvAndExecuteSingleKernel(exeQ, k, 0, nullptr, &k->kernelEvent);
	  k->waitForKernelEvent(exeQ);
	  printMatrix(FLOAT_PTR(bufferC), vDim, vDim);
	  sleep(1);
	  if (i + 1 == aks->getNumKernels()) {
		float sumA = 0;
		float sumC = 0;
		for (int kI = 0; kI < size; ++kI) {
		  sumC += FLOAT_PTR(bufferC)[kI];
		  sumA += FLOAT_PTR(bufferA)[kI];
		}
		flushed_printf("Kernel %s, sliceFactor %d, test:%s ",
					   k->name.c_str(), sliceFactor,
					   sumA == sumC ? "\033[1;32m PASS!\033[0m" : "\033[1;31m FAIL!!!\033[0m");
	  }
	}
	flushed_printf("Requested num chunks: %d, actual: %d\n", pow2(sliceFactor), aks->getNumKernels());
//	delete aks;
	aks = sliceFunc(kernel, ++sliceFactor);
	memset(bufferC->hostPtr, 0, nByte);
//	break;
  }
}

void testSlicingBulk(EDKernel kernel, AtomicKernelSet sliceFunc(EDKernel, uint)) {
  uint sliceFactor = 0;
  auto aks = sliceFunc(kernel, sliceFactor);
  while (aks != nullptr) {
	double timer = getCurrentTime();
	for (int i = 0; i < aks->getNumKernels(); i++) {
	  EDKernel k = aks->kernels[i];
	  edcl->confirmExeEnvAndExecuteSingleKernel(exeQ, k, 0, nullptr, &k->kernelEvent);
//	  k->waitForKernelEvent(exeQ);
	}
	for (int i = 0; i < aks->getNumKernels(); i++) {
	  EDKernel k = aks->kernels[i];
	  k->waitForKernelEvent(exeQ);
	}
	timer = getCurrentTime() - timer;
	float sumA = 0;
	float sumC = 0;
	for (int kI = 0; kI < size; ++kI) {
	  sumC += FLOAT_PTR(bufferC)[kI];
	  sumA += FLOAT_PTR(bufferA)[kI];
	}
	flushed_printf("Kernel %s, sliceFactor %d, test:%s ",
				   kernel->name.c_str(), sliceFactor,
				   sumA == sumC ? "\033[1;32m PASS!\033[0m" : "\033[1;31m FAIL!!!\033[0m");
	flushed_printf("Requested num chunks: %d, actual: %d, time: %f\n",
				   pow2(sliceFactor), aks->getNumKernels(), timer);
//	delete aks;
	aks = sliceFunc(kernel, ++sliceFactor);
	memset(bufferC->hostPtr, 0, nByte);
  }
}

// wrapper functions
AtomicKernelSet edcl_slicingKernelPow2(EDKernel k, uint factor) {
  printf("\033[1;34m%s: ", __func__);
  return edcl->slicingKernelPow2(k, factor);
}

AtomicKernelSet edcl_slicingKernelPow2_extra(EDKernel k, uint factor) {
  printf("\033[1;34m%s: ", __func__);
  return edcl->slicingKernelPow2_extra(k, factor);
}

AtomicKernelSet edcl_slicingKernelPow2_extra_equal(EDKernel k, uint factor) {
  printf("\033[1;34m%s: \n", __func__);
  return edcl->slicingKernelPow2_extra_equal(k, factor);
}
// end wrapper functions

int main() {
  init();
//  testSlicingBulk(vecAdd, slicingKernelPow2);
//  testSlicingBulk(matAdd, slicingKernelPow2);
//  testSlicingBulk(vecAdd, slicingKernelPow2_extra);
//  testSlicingBulk(matAdd, slicingKernelPow2_extra);
//  testSlicingBulk(vecAdd, slicingKernelPow2_extra_equal);
//  testSlicingBulk(matAdd, slicingKernelPow2_extra_equal);
//
//  testSlicingBulk(vecAdd, edcl_slicingKernelPow2);
//  testSlicingBulk(matAdd, edcl_slicingKernelPow2);
//  testSlicingBulk(vecAdd, edcl_slicingKernelPow2_extra);
//  testSlicingBulk(matAdd, edcl_slicingKernelPow2_extra);
  testSlicingBulk(vecAdd, edcl_slicingKernelPow2_extra_equal);
  testSlicingBulk(matAdd, edcl_slicingKernelPow2_extra_equal);

  return 0;
}