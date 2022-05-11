#include "polybench_kernelTest.h"
#include "EDCL_Scheduler.h"

int main() {
  edcl = new EDCL();
  EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
  MM2 mm2(edcl);
  runNative(10, &mm2);
  loopTestKernel(10, scheduler, &mm2);

  MM2 mm2_512(edcl, 512, 512, 512, 512);
  runNative(50, &mm2_512);
  loopTestKernel(50, scheduler, &mm2_512);

  delete edcl;
  return 0;
}
