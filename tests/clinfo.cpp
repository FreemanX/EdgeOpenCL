#include <iostream>
#include "CLRT.h"

int main() {
  printf("GPU: \n");
  auto gpu_cl = CLRT::createGPURunTime();
  gpu_cl->printCLInfo();

  printf("CPU: \n");
  auto pocl_cl = CLRT::createCPURunTime();
  pocl_cl->printCLInfo();

  delete pocl_cl;
  delete gpu_cl;
}
