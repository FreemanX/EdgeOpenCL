#include "EDCL.h"
#include "EDCL_Scheduler.h"

const char *source =
#include "vecAdd.cl"
unsigned int vDim = 4096;
uint n = vDim * vDim;
size_t nByte = n * sizeof(float);
size_t glob_MA[2] = {vDim / 4, vDim};
size_t local_MA[2] = {32, 32};
size_t glob_VA[1] = {n};
size_t local_VA[1] = {1024};

int main() {
  EDCL *edcl = new EDCL();
  // Init Buffers
  EDBuffer A = edcl->createBuffer(nByte);
  EDBuffer B = edcl->createBuffer(nByte);
  EDBuffer C = edcl->createBuffer(nByte);
  randArrayGenerator<float>(0.001, 1, (float *)A->hostPtr, n);
  randArrayGenerator<float>(0.001, 1, (float *)B->hostPtr, n);
  memset(C->hostPtr, 0, nByte);
  // Init Kernels
  EDProgram prog = edcl->createProgram(&source);
  EDKernel matAdd = edcl->createKernel(prog, "matAdd", 5); // program, kernel name, number of args
  EDKernel vecAdd = edcl->createKernel(prog, "vecAdd", 4);
  // Config kernel: work_dim, global_work_size, local_work_size, arg1, arg2...
  matAdd->configKernel(2, glob_MA, local_MA, A, B, C, &vDim, &vDim);
  matAdd->addSubscriber(edcl, vecAdd);
  vecAdd->configKernel(1, glob_VA, local_VA, A, B, C, &n);
  // Get a scheduler instance(Singleton)
  auto scheduler = EDCL_Scheduler::getInstance(edcl);
  scheduler->submitKernels(matAdd, vecAdd);
  // Create scheduling strategy and set strategy parameters
  Sequential strategy(edcl);
  strategy.setExecutionDevice(GPU);
  // Execution
  scheduler->executeKernels(strategy);
  std::cout << "Planing overhead: " << strategy.planingTime;
  std::cout << ", execution time: " << strategy.executionTime << "\n";

  delete edcl;
  return 0;
}