#include <vector>
#include <EDCL.h>
#include "CLRT.h"
#include "polybench.h"
#include "EDCL_Scheduler.h"
#include "polybench_kernelTest.h"

CLRT *cpuRT;
CLRT *gpuRT;
EDQueue cq;
EDQueue gq;
EDCL_Scheduler *scheduler;

void buildProgram(const char **source) {
  cpuRT->createProgramFromSource(source);
  gpuRT->createProgramFromSource(source);
  edcl->createProgram(source);
}

void executeAtomicKernelSet(AtomicKernelSet &kernelSet, EDQueue &queue) {
  for (auto &kernel : kernelSet->kernels) edcl->confirmExeEnv(queue, kernel);
  double timer = getCurrentTime();
  for (auto &kernel : kernelSet->kernels) {
	edcl->executeKernel(queue,
						kernel,
						kernel->num_events_in_wait_list,
						kernel->event_wait_list,
						&kernel->kernelEvent);
  }
  for (auto &kernel : kernelSet->kernels) {
	kernel->waitForKernelEvent(queue);
  }
  timer = getCurrentTime() - timer;
  kernelSet->atomicExeTime = timer;
}

void executePoly(BenchKernel &poly) {
  auto aks = poly.createBenchAKS();
  executeAtomicKernelSet(aks, cq);
  flushed_printf("cq time: %f\n", aks->atomicExeTime);
  executeAtomicKernelSet(aks, gq);
  flushed_printf("gq time: %f\n", aks->atomicExeTime);
}

void test2mm() {
  const char *source_2mm =
#include "Kernels/2mm.cl"
  buildProgram(&source_2mm);
  MM2 mm_2(edcl, 128, 128, 128, 128);
  executePoly(mm_2);

  //scheduler->clearKernels();
  //loopTestKernel(10, scheduler, &mm_2);
}

int main() {
  cpuRT = CLRT::createCPURunTime();
  gpuRT = CLRT::createGPURunTime();
  edcl = new EDCL();
  cq = edcl->createDeviceCmdQueueProfilingEnabled(CPUQueue);
  gq = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  scheduler = GET_SCHEDULER(edcl); // TODO: runtime problem
  test2mm();

  return 0;
}