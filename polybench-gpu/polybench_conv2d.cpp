#include "polybench_kernelTest.h"
#include "EDCL_Scheduler.h"

int main() {
  edcl = new EDCL();
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);

  Conv2D conv2d(edcl);
  runNative(100, &conv2d);
  loopTestKernel(100, scheduler, &conv2d);

  Conv2D conv2d_1024(edcl, 1024, 1024);
  runNative(100, &conv2d_1024);
  loopTestKernel(100, scheduler, &conv2d_1024);

  delete edcl;
  return 0;
}
