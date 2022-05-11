#include "polybench_stgTest.h"
#include "Stg_CA2020.h"

int main() {
//  INIT_TEST;
  INIT_TEST

  CA2020 ca2020(edcl);
  scheduler->executeKernels(ca2020);
  double overallTime;
  overallTime = scheduler->executeKernels(ca2020);
  flushed_printf("\tOverall time: %lf\n", overallTime);
  printExeTime(scheduler);
  return 0;
}

