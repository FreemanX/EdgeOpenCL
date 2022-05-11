#include "polybench_kernelTest.h"
#include "EDCL_Scheduler.h"

int main() {
  edcl = new EDCL();
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);

  Conv3D conv3d(edcl);
  runNative(10, &conv3d);
  loopTestKernel(10, scheduler, &conv3d);

  Conv3D conv3d_128(edcl, 128, 128, 128);
  runNative(10, &conv3d_128);
  loopTestKernel(10, scheduler, &conv3d_128);

  delete edcl;
  return 0;
}
