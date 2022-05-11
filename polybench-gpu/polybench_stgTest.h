#ifndef EDCL_POLYBENCH_GPU_POLYBENCH_STGTEST_H_
#define EDCL_POLYBENCH_GPU_POLYBENCH_STGTEST_H_
#include "CPU_utils.h"
#include "polybench.h"
#include "EDCL_Scheduler.h"

#define INIT_TEST                                                              \
  EDCL *edcl = new EDCL();                                                     \
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);                             \
  Conv2D conv2d(edcl);                                                         \
  Conv3D conv3d(edcl);                                                         \
  MM2 mm2(edcl);                                                               \
  ATAX atax(edcl);                                                             \
  BICG bicg(edcl);                                                             \
  Correlation corr(edcl);                                                      \
  Covariance covar(edcl);                                                      \
  FDTD2D fdtd2d(edcl);                                                         \
  GEMM gemm(edcl);                                                             \
  Gesummv gesummv(edcl);                                                       \
  GramSchmidt gram_schmidt(edcl);                                              \
  mvt mvt(edcl);                                                               \
  syr syrk(edcl);                                                              \
  scheduler->submitAtomicKernelSets(                                           \
      conv2d.createBenchAKS(), conv3d.createBenchAKS(), mm2.createBenchAKS(),  \
      atax.createBenchAKS(), bicg.createBenchAKS(), corr.createBenchAKS(),     \
      covar.createBenchAKS(), fdtd2d.createBenchAKS(), gemm.createBenchAKS(),  \
      gesummv.createBenchAKS(), gram_schmidt.createBenchAKS(),                 \
      mvt.createBenchAKS(), syrk.createBenchAKS());

#define RUN_NATIVE                                                             \
  {                                                                            \
    double nativeOverallTime;                                                  \
    double timer[13];                                                          \
    int numThreads[4] = {1, 2, 4, 8};                                          \
    for (int n : numThreads) {                                                 \
      if (n == 1) while (setThreadAffinity(7) != 0) {}               \
      flushed_printf("Running polybench native, num threads = %d...\n", n);    \
      nativeOverallTime = getCurrentTime();                                    \
      timer[0] = conv2d.nativeKernel(n);                                       \
      timer[1] = conv3d.nativeKernel(n);                                       \
      timer[2] = mm2.nativeKernel(n);                                          \
      timer[3] = atax.nativeKernel(n);                                         \
      timer[4] = bicg.nativeKernel(n);                                         \
      timer[5] = corr.nativeKernel(n);                                         \
      timer[6] = covar.nativeKernel(n);                                        \
      timer[7] = fdtd2d.nativeKernel(n);                                       \
      timer[8] = gemm.nativeKernel(n);                                         \
      timer[9] = gesummv.nativeKernel(n);                                      \
      timer[10] = gram_schmidt.nativeKernel(n);                                \
      timer[11] = mvt.nativeKernel(n);                                         \
      timer[12] = syrk.nativeKernel(n);                                        \
      nativeOverallTime = getCurrentTime() - nativeOverallTime;                \
      flushed_printf("\tOverall time: %lf\n", nativeOverallTime);              \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     conv2d.getName(), timer[0]);                              \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     conv3d.getName(), timer[1]);                              \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     mm2.getName(), timer[2]);                                 \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     atax.getName(), timer[3]);                                \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     bicg.getName(), timer[4]);                                \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     corr.getName(), timer[5]);                                \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     covar.getName(), timer[6]);                               \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     fdtd2d.getName(), timer[7]);                              \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     gemm.getName(), timer[8]);                                \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     gesummv.getName(), timer[9]);                             \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     gram_schmidt.getName(), timer[10]);                       \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     mvt.getName(), timer[11]);                                \
      flushed_printf("\t(numThreads=%d) %s native time: %lf\n", n,             \
                     syrk.getName(), timer[12]);                               \
    }                                                                          \
  }

EXECUTE_DEVICE devices[8] = {GPU, CPU, B_SD4, B_SD2, B_SD1, L_SD4, L_SD2, L_SD1};

void printExeTime(EDCL_Scheduler *scheduler) {
	for (const auto &ks : scheduler->atomicKernelSets) {
		flushed_printf("\t  AKS %d time = %lf s\n", ks->getTrackNum(),
									 ks->atomicExeTime);
		// if number of kernels in a kernel set greater than 10, print each kernel
		// time otherwise, print average kernel time
		if (ks->getNumKernels() < 10) {
			for (auto &k : ks->kernels) {
				flushed_printf("\t    Kernel %s, clEvent time = %lf\n", k->name.c_str(),
											 k->getKernelEventTime());
			}
		} else {
			double totalTime = 0;
			for (auto &k : ks->kernels) {
				totalTime += k->getKernelEventTime();
			}
			flushed_printf("\t    Kernel average clEvent time = %lf\n",
										 totalTime / ks->getNumKernels());
		}
	}
}

#endif // EDCL_POLYBENCH_GPU_POLYBENCH_STGTEST_H_
