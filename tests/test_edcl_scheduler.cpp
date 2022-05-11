#include "Stg_Sequential.h"
#include "EDCL_Scheduler.h"
#include "EDCL.h"
#include "testlib_cpu_vector_add.h"
#include "testlib_macros.h"

#define FLOAT_PTR(EDBuffer) HOST_PTR(float, EDBuffer)
const char *source =
#include "vecAdd.cl"
unsigned int vDim = 4096;
uint size = vDim * vDim;
size_t nByte = size * sizeof(float);
double timer;

EXECUTE_DEVICE devices[8] = {L_SD1, L_SD2, L_SD4, B_SD1, B_SD2, B_SD4, GPU, CPU};

void printAvgExeTime(EDCL_Scheduler *scheduler) {
	double totalKernelTime = 0;
	double totalAKSTime = 0;
	for (const auto &ks : scheduler->atomicKernelSets) {
		totalAKSTime += ks->atomicExeTime;
		for (auto k : ks->kernels)
			totalKernelTime += k->getKernelEventTime();
	}
	flushed_printf("\tKernel avg execution time: %lf\n",
								 totalKernelTime / scheduler->getNumEDKernels());
}

int main() {
	EDCL *edcl = new EDCL();
	int loopNum = 200;
	// Init buffers
	EDBuffer bufferA = edcl->createBuffer(nByte);
	EDBuffer bufferB = edcl->createBuffer(nByte);
	TIME_CODE_ONCE(
			timer,
			randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferA), size);
					randArrayGenerator<float>(0.001, 1, FLOAT_PTR(bufferB), size);
	)
	EDBuffer bufferC = edcl->createBuffer(nByte);
	memset(bufferC->hostPtr, 0, nByte);

	auto *c_verify = static_cast<float *>(malloc(nByte));
	memset(c_verify, 0, nByte);
	vecAddOptimized(FLOAT_PTR(bufferA), FLOAT_PTR(bufferB), c_verify, size);

	// Init Kernels
	EDProgram program = edcl->createProgram(&source);
	size_t vecAdd_globalSize[1] = {size};
	size_t matAdd_globalSize[2] = {vDim / 4, vDim};

	// Execute with scheduler
	auto scheduler = GET_SCHEDULER(edcl);
	Sequential sequential(edcl);
	double overallTime;

//  for (int i = 0; i < loopNum; ++i) {
//	EDKernel vecAdd = edcl->createKernel(program, "vecAdd", 4);
//	vecAdd->configKernel(1, vecAdd_globalSize, nullptr, bufferA, bufferB, bufferC, &size);
//	scheduler->submitKernels(vecAdd);
//  }
//  for (auto &device : devices) { // test every device
//	flushed_printf("Running vecAdd on %s\n", EXECUTE_DEVICE_ToString(device));
//	sequential.setExecutionDevice(device);
//	overallTime = scheduler->executeKernels(sequential);
//	flushed_printf("\tOverall time: %lf\n", overallTime);
//	printAvgExeTime(scheduler);
//  }
//
//  scheduler->clearKernels();
//
//  for (int i = 0; i < loopNum; ++i) {
//	EDKernel matAdd = edcl->createKernel(program, "matAdd", 5);
//	matAdd->configKernel(2, matAdd_globalSize, nullptr, bufferA, bufferB, bufferC, &vDim, &vDim);
//	scheduler->submitKernels(matAdd);
//  }
//  for (auto &device : devices) { // test every device
//	flushed_printf("Running matAdd on %s\n", EXECUTE_DEVICE_ToString(device));
//	sequential.setExecutionDevice(device);
//	overallTime = scheduler->executeKernels(sequential);
//	flushed_printf("\tOverall time: %lf\n", overallTime);
//	printAvgExeTime(scheduler);
//  }

	scheduler->clearKernels();

	for (int i = 0; i < loopNum; ++i) {
		EDKernel matAdd = edcl->createKernel(program, "matAdd", 5);
		matAdd->configKernel(2, matAdd_globalSize, nullptr, bufferA, bufferB, bufferC, &vDim, &vDim);
		scheduler->submitKernels(matAdd);
	}

	flushed_printf("\nRunning matAdd on %s using individual AKS\n", "GPU");
	sequential.setExecutionDevice(GPU);
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);
	flushed_printf("Running matAdd on %s using individual AKS\n", "CPU");
	sequential.setExecutionDevice(B_SD4);
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);

	scheduler->clearKernels();

	AtomicKernelSet vecAddSets = edcl->createAtomicKernelSet(loopNum);
	for (int i = 0; i < loopNum; ++i) {
		EDKernel matAdd = edcl->createKernel(program, "matAdd", 5);
		matAdd->configKernel(2, matAdd_globalSize, nullptr, bufferA, bufferB, bufferC, &vDim, &vDim);
		vecAddSets->addKernel(matAdd);
	}
	scheduler->submitAtomicKernelSets(vecAddSets);
	flushed_printf("\nRunning matAdd on %s using one AKS\n", "GPU");
	sequential.setExecutionDevice(GPU);
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);
	flushed_printf("Running matAdd on %s using one AKS\n", "CPU");
	sequential.setExecutionDevice(B_SD4);
	overallTime = scheduler->executeKernels(sequential);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);

	scheduler->clearKernels();

	SequentialImproved sequential_improved(edcl);

	for (int i = 0; i < loopNum; ++i) {
		EDKernel matAdd = edcl->createKernel(program, "matAdd", 5);
		matAdd->configKernel(2, matAdd_globalSize, nullptr, bufferA, bufferB, bufferC, &vDim, &vDim);
		scheduler->submitKernels(matAdd);
	}

	flushed_printf("\nRunning matAdd on %s using individual AKS under improved sequential strategy\n", "GPU");
	sequential_improved.setExecutionDevice(GPU);
	overallTime = scheduler->executeKernels(sequential_improved);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);
	flushed_printf("Running matAdd on %s using individual AKS under improved sequential strategy\n", "CPU");
	sequential_improved.setExecutionDevice(B_SD4);
	overallTime = scheduler->executeKernels(sequential_improved);
	flushed_printf("\tOverall time: %lf\n", overallTime);
	printAvgExeTime(scheduler);

	delete edcl;
	return 0;
}
