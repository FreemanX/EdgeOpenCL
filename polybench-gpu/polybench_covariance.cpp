#include "polybench_kernelTest.h"
#include "EDCL_Scheduler.h"

#define BENCH_KERNEL Covariance

int main() {
  edcl = new EDCL();
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);

  BENCH_KERNEL defaultBench(edcl);
  runNative(10, &defaultBench);
  loopTestKernel(10, scheduler, &defaultBench);

  BENCH_KERNEL changedBench(edcl, 512, 512);
  loopTestKernel(10, scheduler, &changedBench);
  runNative(10, &changedBench);

  delete edcl;
  return 0;
}



