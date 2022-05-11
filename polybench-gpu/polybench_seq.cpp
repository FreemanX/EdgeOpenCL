#include <unistd.h>
#include "polybench.h"
#include "EDCL_Scheduler.h"
#include "polybench_stgTest.h"

int main() {
  INIT_TEST
  RUN_NATIVE

// Strategy test
  Sequential sequential(edcl);
  // Dry run
  sequential.setExecutionDevice(CPU);
  scheduler->executeKernels(sequential);
  // Dry run
  sequential.setExecutionDevice(GPU);
  scheduler->executeKernels(sequential);

  double overallTime;
  for (auto &device : devices) { // test on every device
	sleep(60);
	sequential.setExecutionDevice(device);
	flushed_printf("Running polybench on %s\n",
				   EXECUTE_DEVICE_ToString(device));
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printExeTime(scheduler);
  }

  delete edcl;
  return 0;
}
