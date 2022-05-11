#ifndef EDCL_POLYBENCH_GPU_POLYBENCH_KERNELTEST_H_
#define EDCL_POLYBENCH_GPU_POLYBENCH_KERNELTEST_H_

#include <EDCL_Scheduler.h>
#include "polybench.h"

EDCL *edcl;

//EXECUTE_DEVICE devices[8] = {L_SD4, L_SD1, L_SD2, B_SD1, B_SD2, B_SD4, GPU, CPU};
EXECUTE_DEVICE devices[8] = {GPU, CPU, B_SD4, B_SD2, B_SD1, L_SD4, L_SD2, L_SD1,};
//EXECUTE_DEVICE devices[4] = {GPU, CPU, B_SD4, L_SD4};
//EXECUTE_DEVICE devices[2] = {CPU, GPU};
void printAvgExeTime(EDCL_Scheduler *scheduler) {
  double totalAKSTime = 0;
  int numAKS = scheduler->atomicKernelSets.size();
  int numKernels = scheduler->atomicKernelSets[0]->getNumKernels(); // num kernels in each AKS
  if (numKernels > 5) numKernels = 5;
  double totalKernelTime[numKernels];
  for (const auto& ks : scheduler->atomicKernelSets) {
	totalAKSTime += ks->atomicExeTime;
	for (int i = 0; i < numKernels; ++i) {
	  totalKernelTime[i] += ks->kernels[i]->getKernelEventTime();
	}
  }
  flushed_printf("\tAKS avg time: %lf \n", totalAKSTime / numAKS);
  for (int i = 0; i < numKernels; i++) {
	flushed_printf("\tKernel %s avg execution time: %lf\n",
				   scheduler->atomicKernelSets[0]->kernels[i]->name.c_str(),
				   totalKernelTime[i] / numAKS);
  }
}

void loopTestKernel(int loopNum, EDCL_Scheduler *scheduler, BenchKernel *kernel) {
  Sequential sequential(edcl);
  for (int i = 0; i < loopNum; ++i)
	scheduler->submitAtomicKernelSets(kernel->createBenchAKS());
  double overallTime;
  for (auto &device : devices) { // test on every device
	flushed_printf("Running %s on %s\n", kernel->getName(), EXECUTE_DEVICE_ToString(device));
	sequential.setExecutionDevice(device);
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);
  }
  scheduler->clearKernels();
}

void runNative(int loopNum, BenchKernel *kernel) {
  int numThreads[4] = {1, 2, 4, 8};
  for (int n : numThreads) {
	double totalTime = 0;
	for (int i = 0; i < loopNum; ++i) {
	  totalTime += kernel->nativeKernel(n);
	}
	flushed_printf("===Native %s (numThreads=%d) avg time: %lf===\n", kernel->getName(), n, totalTime / loopNum);
  }
}
#endif //EDCL_POLYBENCH_GPU_POLYBENCH_KERNELTEST_H_
