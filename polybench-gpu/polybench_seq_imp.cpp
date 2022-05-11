#include <unistd.h>
#include "polybench.h"
#include "EDCL_Scheduler.h"
#include "polybench_stgTest.h"

int main() {
  INIT_TEST
  RUN_NATIVE

  // Strategy test
  SequentialImproved sequentialImproved(edcl);
  // Dry run
  sequentialImproved.setExecutionDevice(CPU);
  scheduler->executeKernels(sequentialImproved);
  // Dry run
  sequentialImproved.setExecutionDevice(GPU);
  scheduler->executeKernels(sequentialImproved);

  double overallTime;
  for (auto &device : devices) { // test on every device
	sleep(60);
	sequentialImproved.setExecutionDevice(device);
	flushed_printf("Running polybench on %s\n",
				   EXECUTE_DEVICE_ToString(device));
	overallTime = scheduler->executeKernels(sequentialImproved);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printExeTime(scheduler);
  }

  delete edcl;
  return 0;
}
