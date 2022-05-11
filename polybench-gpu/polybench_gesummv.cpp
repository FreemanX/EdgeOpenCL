#include "polybench_kernelTest.h"
#include "EDCL_Scheduler.h"

#define BENCH_KERNEL Gesummv

int main() {
  edcl = new EDCL();
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);

  BENCH_KERNEL defaultBench(edcl);
  runNative(10, &defaultBench);
  loopTestKernel(10, scheduler, &defaultBench);

  BENCH_KERNEL changedBench(edcl, 1024);
  runNative(10, &changedBench);
  loopTestKernel(10, scheduler, &changedBench);

  delete edcl;
  return 0;
}

