#include <tar.h>
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
/*
 * Deprecated! Actual implementation using test_clrt_each_subdevice.cpp
 * */
#include <vector>
#include "CLRT.h"
#include "utils.h"
#include "EDCL.h"
#include "testlib_macros.h"
#include "testlib_cpu_vector_add.h"
#include "CPU_utils.h"
#include "pthread.h"
#define FLOAT_PTR(EDBuffer) HOST_PTR(float, EDBuffer)

// Test functions
uint numSubDeviceDivisions = 0;
// Array of sub-device vectors. Each vector contains sub-devices with the same number of CUs.
std::vector<std::vector<cl_device_id> *> sub_device_collection; // CLRT
std::vector<std::vector<cl_command_queue> *> sub_device_queue_collection; // CLRT
void initSubDevicesAndQueues(CLRT *poclRT); // CLRT
// NumCU can be 1, 2, 4, 8...
EDQueue getSubDeviceExeQueue(CLRT *clrt, uint NumCU); // EDCL
void releaseSubDeviceExeQueue(EDQueue Queue); // EDCL
// Test functions


// Test
unsigned int vDim = 64;
uint size = vDim * vDim;
size_t nByte = size * sizeof(float);
const char *source =
#include "vecAdd.cl"
EDBuffer vecA;
EDBuffer vecB;
EDBuffer vecC;
EDBuffer matC;
float *c_verify;
EDProgram program;
EDKernel vecAdd;
EDKernel matAdd;
size_t vecAdd_globalSize[] = {size};
size_t matAdd_globalSize[] = {vDim / 4, vDim};
void initTest(EDCL *edcl);
void endTest(EDCL *edcl);
void multiSubDeviceQueueTest(EDCL *edcl);
void *threadWorkFunction(void *WP);// WP --> work pack

int main() {
  EDCL *edcl = new EDCL();
  initSubDevicesAndQueues(edcl->cpuRT_);
  initTest(edcl);
  multiSubDeviceQueueTest(edcl);
  endTest(edcl);
  delete edcl;
  return 0;
}

typedef struct InternalPThreadWorkPack_ {
  EDCL *edcl;
  EDKernel kernel;
  EDQueue queue;
  __unused InternalPThreadWorkPack_(EDCL *EDCL, EDKernel Kernel, EDQueue Queue) {
	edcl = EDCL;
	kernel = Kernel;
	queue = Queue;
  }
} InternalPThreadWorkPack_;

typedef InternalPThreadWorkPack_ *WorkPack;

void multiSubDeviceQueueTest(EDCL *edcl) {
//  cl_int err;
//  setThreadAffinity(0);
//  setThreadAffinity({4, 5, 6, 7});
  uint numCUs[] = {1, 2, 4};
//  vecAddOMP8Opt(FLOAT_PTR(vecA), FLOAT_PTR(vecB), c_verify, 8);
  for (auto numCU : numCUs) {
	flushed_printf("Num CU: %d\n", numCU);
	EDQueue vecAddQ = getSubDeviceExeQueue(edcl->cpuRT_, numCU);
	EDQueue matAddQ = getSubDeviceExeQueue(edcl->cpuRT_, numCU);
	int cpu_id1, cpu_id2;
	for (int i = 0; i < numCU; ++i) {
	  cpu_id1 = 7 - i;
	  cpu_id2 = 7 - (int)numCU - i;
	  SET_BIT(vecAddQ->cpuSet, cpu_id2);
	  SET_BIT(matAddQ->cpuSet, cpu_id1);
	}
//	cl_event syncEvent = clCreateUserEvent(edcl->cpuRT_->context_, &err);
//	CLRT_ERR(err);
	auto vecAddPack = new InternalPThreadWorkPack_(edcl, vecAdd, vecAddQ);
	auto matAddPack = new InternalPThreadWorkPack_(edcl, matAdd, matAddQ);
	pthread_t vecAddThread, matAddThread;
	double *exeTime[2];
	pthread_create(&vecAddThread, nullptr, threadWorkFunction, vecAddPack);
	pthread_create(&matAddThread, nullptr, threadWorkFunction, matAddPack);
	pthread_join(vecAddThread, (void **)&exeTime[0]);
	pthread_join(matAddThread, (void **)&exeTime[1]);
	flushed_printf("  vecAdd exeTime: %lf, matAdd exeTime: %lf\n",
				   *exeTime[0], *exeTime[1]);
	flushed_printf("  vecAdd memcmp: %d, matAdd memcmp: %d\n",
				   memcmp(c_verify, vecC->hostPtr, nByte), memcmp(c_verify, matC->hostPtr, nByte));
	releaseSubDeviceExeQueue(vecAddQ);
	releaseSubDeviceExeQueue(matAddQ);
	free(exeTime[0]);
	free(exeTime[1]);
	delete vecAddPack;
	delete matAddPack;
	memset(vecC->hostPtr, 0, nByte);
	memset(matC->hostPtr, 0, nByte);
  }
}

void *threadWorkFunction(void *WP) {
  // Get work pack
  auto wp = (WorkPack)WP;
  auto *timer = static_cast<double *>(calloc(1, sizeof(double)));
  // set CPU affinity
  while (setThreadAffinity(wp->queue->cpuSet) != 0) {}
// dry run
  cl_event event;
  wp->edcl->confirmExeEnvAndExecuteSingleKernel(wp->queue, wp->kernel, 0, nullptr, &event);
  clWaitForEvents(1, &event, wp->queue->clrt->oclLibHandle_);
  wp->queue->finish();
  // set CPU affinity
  while (setThreadAffinity(wp->queue->cpuSet) != 0) {}
  printThreadAffinity(wp->kernel->numArgs);
  // execute kernel
  TIME_CODE_LOOP(*timer, 10,
				 wp->edcl->executeKernel(wp->queue, wp->kernel, 0, nullptr, &event);
					 clWaitForEvents(1, &event, wp->queue->getSoHandle());
  )
  wp->queue->finish();

  flushed_printf("Kernel %s, event timer: %lf\n",
				 wp->kernel->name.c_str(), wp->queue->clrt->getExeTime(&event));
  // exit
  pthread_exit(timer);
}

void endTest(EDCL *edcl) {
  edcl->releaseEDProgram(program);
  delete vecAdd;
  edcl->releaseEDBuffer(vecA);
  edcl->releaseEDBuffer(vecB);
  edcl->releaseEDBuffer(vecC);
  edcl->releaseEDBuffer(matC);
  delete c_verify;
}

void initTest(EDCL *edcl) {
  double timer;
  vecA = edcl->createBuffer(nByte);
  vecB = edcl->createBuffer(nByte);
  flushed_printf("Initializing vector(%d)... ", size);
  TIME_CODE_ONCE(
	  timer,
	  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(vecA), size);
		  randArrayGenerator<float>(0.001, 1, FLOAT_PTR(vecB), size);
  )
  flushed_printf("done. Time cost: %lf\n", timer);
  vecC = edcl->createBuffer(nByte);
  matC = edcl->createBuffer(nByte);
  c_verify = static_cast<float *>(malloc(nByte));
  memset(vecC->hostPtr, 0, nByte);
  memset(matC->hostPtr, 0, nByte);
  memset(c_verify, 0, nByte);
  TIME_CODE_ONCE(timer,
				 vecAddOptimized(FLOAT_PTR(vecA), FLOAT_PTR(vecB), c_verify, size);)
  flushed_printf("Ground truth calculated, timer: %lf\n", timer);
  program = edcl->createProgram(&source);
  vecAdd = edcl->createKernel(program, "vecAdd", 4);
  matAdd = edcl->createKernel(program, "matAdd", 5);
  vecAdd->configKernel(1, vecAdd_globalSize, nullptr, vecA, vecB, vecC, &size);
  matAdd->configKernel(2, matAdd_globalSize, nullptr, vecA, vecB, matC, &vDim, &vDim);
}

// EDCL func
EDQueue getSubDeviceExeQueue(CLRT *CLRT, uint NumCU) {
  auto q = new InternalEDQueue_(CLRT);
  q->numCU = NumCU;
  q->collectionIdx = log2(NumCU);
  q->queueType = CPUQueue;
  q->queue = sub_device_queue_collection[q->collectionIdx]->back();
  sub_device_queue_collection[q->collectionIdx]->pop_back();
  return q;
}

// EDCL func
void releaseSubDeviceExeQueue(EDQueue Queue) {
  sub_device_queue_collection[Queue->collectionIdx]->push_back(Queue->queue);
  delete Queue;
}

// CLRT cpu instance func
void initSubDevicesAndQueues(CLRT *poclRT) {
  uint numCUPerDevice = poclRT->getMaxSubDevice(0, 0);
  std::cout << "Init numCUPerDevice " << numCUPerDevice << "\n";
  while (numCUPerDevice > 0) {
	numSubDeviceDivisions++;
	numCUPerDevice = numCUPerDevice >> 1;
	std::cout << "numCUPerDevice " << numCUPerDevice << "\n";
  }
  numCUPerDevice = 1; // reset to 1
  std::cout << "numSubDeviceDivisions " << numSubDeviceDivisions << "\n";
  uint *subDeviceDivisions = (uint *)calloc(numSubDeviceDivisions, sizeof(uint));
  for (int kI = 0; kI < numSubDeviceDivisions; ++kI) {
	subDeviceDivisions[kI] = numCUPerDevice << kI;
	std::cout << "subDeviceDivisions[" << kI << "] " << subDeviceDivisions[kI] << "\n";
  }
  cl_device_partition_property props[3];
  props[0] = CL_DEVICE_PARTITION_EQUALLY;  // Equally
  props[2] = 0;                            // End of the property list

  for (int i = 0; i < numSubDeviceDivisions; ++i) {
	numCUPerDevice = subDeviceDivisions[i];
	props[1] = numCUPerDevice;     // compute units per sub-device
	// Specifies the size of the out_devices array: max num of sub-devices(num cpus) / num of cu(cpu core) per device
	uint numSubDevices = poclRT->getMaxSubDevice(0, 0) / numCUPerDevice;
// Provides a buffer for the generated subdevices with a number
// of elements specified by num_sub_devices
//	auto *devices = new std::vector<cl_device_id>(numSubDevices);
	sub_device_collection.push_back(new std::vector<cl_device_id>(numSubDevices));
	sub_device_queue_collection.push_back(new std::vector<cl_command_queue>);
	std::cout << "Size of sub_device_collection: " << sub_device_collection.size() << "\n";
	std::cout << "Size of sub_device_queue_collection: " << sub_device_queue_collection.size() << "\n";
// clCreateSubDevices returns the number of subdevices
// in which the device may be partitioned into considering the
// partition type and the other values specified in the property list
	uint numSubDeviceRet = 0;
	CLRT_ERR(clCreateSubDevices(poclRT->device_id_,
								props,
								numSubDevices,
								sub_device_collection[i]->data(),
								&numSubDeviceRet, poclRT->oclLibHandle_));
	std::cout << "Num Sub-devices: " << numSubDevices << ", num sub-devices returned: " << numSubDeviceRet << "\n";
	std::cout << "Size of sub_device_collection[" << i << "]: " << sub_device_collection[i]->size() << "\n";
	uint minSubDev = MIN(numSubDevices, numSubDeviceRet);
	std::cout << minSubDev << " -- minSubDev\n";
	for (int j = 0; j < minSubDev; ++j) {
	  std::cout << " Creating queue" << j << "(" << numCUPerDevice << ").\n";
	  cl_command_queue queue = poclRT->
		  createCLExecutionQueue(sub_device_collection[i]->at(j), CL_QUEUE_PROFILING_ENABLE)->command_queue;
	  sub_device_queue_collection[i]->push_back(queue);
//	  clReleaseCommandQueue(queue);
	}
  }
  delete subDeviceDivisions;
  std::cout << "-----------\n";
  std::cout << "| Summery |\n";
  std::cout << "-----------\n";
  std::cout << " numSubDeviceDivisions: \t\t " << numSubDeviceDivisions << "\n";
  std::cout << " Size of Sub-device_collection: \t " << sub_device_collection.size() << "\n";
  std::cout << " Size of Sub-device_queue_collection: \t " << sub_device_queue_collection.size() << "\n";
  for (int i = 0; i < numSubDeviceDivisions; ++i) {
	std::cout << "  Size of sub-device_collection[" << i << "]\t " <<
			  sub_device_collection[i]->size() << "\n";
	std::cout << "  Size of sub-device_queue_collection[" << i << "] " <<
			  sub_device_queue_collection[i]->size() << "\n";
  }
}
