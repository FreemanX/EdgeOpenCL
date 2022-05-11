#include <unistd.h>
#include "EDCL_ResourceManager.h"
#include "EDCL.h"
#include "utils.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
EDCL *edcl;
ResourceManager rm;

#define CPU_STATUS rm->getAVAILABLE_CPU()
#define numCPU rm->getNumCPU()

void bindEDQueueCPU(EDQueue queue) {
  flushed_printf("Binding queue(%d), hetero level %d with %d cores...\n",
				 queue->coreStart, rm->getEDQueueHeteroLevel(queue), queue->numCU);
  flushed_printf("\tCPU_STATUS before bind:\t\t");
  printBits(CPU_STATUS, numCPU);
  flushed_printf("\tqueue->numCU: %d, queue->CPUSet:\t", queue->numCU);
  printBits(queue->cpuSet, numCPU);

  rm->bindEDQueueCPU(queue);

  flushed_printf("\tCPU_STATUS after bind:\t\t");
  printBits(CPU_STATUS, numCPU);
}

void unbindEDQueueCPU(EDQueue queue) {
  flushed_printf("Unbinding queue(%d) with %d cores...\n", queue->coreStart, queue->numCU);
  flushed_printf("\tCPU_STATUS before unbind:\t");
  printBits(CPU_STATUS, numCPU);
  flushed_printf("\tqueue->numCU: %d, queue->CPUSet:\t", queue->numCU);
  printBits(queue->cpuSet, numCPU);

  rm->unbindEDQueueCPU(queue);

  flushed_printf("\tCPU_STATUS after unbind:\t");
  printBits(CPU_STATUS, numCPU);
}

void coreBindingTests() {
  // Core binding tests
  EDQueue b_q1 = edcl->createSubDeviceExeQueue(0, 7);
  EDQueue b_q2 = edcl->createSubDeviceExeQueue(1, 2);
  EDQueue b_q3 = edcl->createSubDeviceExeQueue(2, 1);
  EDQueue l_q1 = edcl->createSubDeviceExeQueue(0, 3);
  EDQueue l_q2 = edcl->createSubDeviceExeQueue(1, 1);
  EDQueue l_q3 = edcl->createSubDeviceExeQueue(2, 0);
  flushed_printf("Test single queue binding and unbinding...\n");
  bindEDQueueCPU(b_q1);
  unbindEDQueueCPU(b_q1);
  bindEDQueueCPU(b_q2);
  unbindEDQueueCPU(b_q2);
  bindEDQueueCPU(b_q3);
  unbindEDQueueCPU(b_q3);
  bindEDQueueCPU(l_q1);
  unbindEDQueueCPU(l_q1);
  bindEDQueueCPU(l_q2);
  unbindEDQueueCPU(l_q2);
  bindEDQueueCPU(l_q3);
  unbindEDQueueCPU(l_q3);

  flushed_printf("Test multiple queue binding and unbinding...\n");
  flushed_printf("Core binding conflict\n");
  bindEDQueueCPU(b_q1);
  bindEDQueueCPU(b_q2);
  bindEDQueueCPU(b_q3);
  bindEDQueueCPU(l_q1);
  bindEDQueueCPU(l_q2);
  bindEDQueueCPU(l_q3);
  unbindEDQueueCPU(b_q1);
  unbindEDQueueCPU(b_q2);
  unbindEDQueueCPU(l_q1);
  unbindEDQueueCPU(l_q2);
  bindEDQueueCPU(b_q3);
  bindEDQueueCPU(l_q3);
  bindEDQueueCPU(b_q3);
  bindEDQueueCPU(l_q3);
  unbindEDQueueCPU(b_q3);
  unbindEDQueueCPU(l_q3);
  unbindEDQueueCPU(b_q1);
  unbindEDQueueCPU(b_q2);
  unbindEDQueueCPU(b_q3);
  unbindEDQueueCPU(l_q1);
  unbindEDQueueCPU(l_q2);
  unbindEDQueueCPU(l_q3);
  edcl->releaseEDQueue(b_q1);
  edcl->releaseEDQueue(b_q2);
  edcl->releaseEDQueue(b_q3);
  edcl->releaseEDQueue(l_q1);
  edcl->releaseEDQueue(l_q2);
  edcl->releaseEDQueue(l_q3);
}

void printResourceInfo() {
  // Basic info
  flushed_printf("Total CPU cores: %d, HeteroLevel: %d. Details:\n",
				 rm->getNumCPU(), rm->getCPUHeteroLevel());
  std::vector<uint> coreCnts = rm->getHeteroCoreNums();
  for (int kI = 0; kI < coreCnts.size(); ++kI)
	flushed_printf("\tHeteroLevel(%d): %d\n", kI, coreCnts[kI]);
  flushed_printf("\tnumBigCores: %d, numLittleCores: %d\n",
				 rm->getNumBigCores(), rm->getNumLittleCores());
  auto cpuFreqs = rm->getHeteroCPUFrequencies();
  for (auto &h_level : cpuFreqs) {
	flushed_printf("\tCPU available frequencies(%d): ", h_level.size());
	for (auto &freq : h_level) flushed_printf("%d ", freq);
	flushed_printf("\n");
  }

  auto gpuFreqs = rm->getGPUFrequencies();
  flushed_printf("\tGPU available frequencies(%d): ", gpuFreqs.size());
  for (auto &freq : gpuFreqs) flushed_printf("%d ", freq);
  flushed_printf("\n");
  rm->setSoCMaxFrequency();
  for (int kJ = 0; kJ < rm->getNumCPU(); ++kJ)
	flushed_printf("\tCPU %d freq: %d\n", kJ, rm->getCPUCurrentFrequency(kJ));
  flushed_printf("\tGPU freq: %d\n", rm->getGPUCurrentFrequency());
}

void testRequestBigCoreQueue() {
  uint numQ = 5;
  EDQueue bigQueues[numQ];
  flushed_printf("Testing creation of bigQ...\n");

  for (int Log2NumCU = 0; Log2NumCU < 4; ++Log2NumCU) {
	flushed_printf("\tInit CPU_STATUS:\t");
	printBits(CPU_STATUS, numCPU);
	for (int kI = 0; kI < numQ; ++kI) {
	  bigQueues[kI] = rm->createUniqQueueBigCore(Log2NumCU, edcl);
	  flushed_printf("\tCPU_STATUS after create bigQ(%d):\t", (int)exp2(Log2NumCU));
	  printBits(CPU_STATUS, numCPU);
	  if (bigQueues[kI] == nullptr) break;
	}

	for (int kI = 0; kI < numQ; ++kI) {
	  if (bigQueues[kI] != nullptr) {
		rm->unbindEDQueueCPU(bigQueues[kI]);
		edcl->releaseEDQueue(bigQueues[kI]);
		bigQueues[kI] = nullptr;
	  } else {
		flushed_printf("\tbigQ(%d) ends at %d\n", (int)exp2(Log2NumCU), kI);
		break;
	  }
	}
  }
}

void testRequestLittleCoreQueue() {
  uint numQ = 5;
  EDQueue bigQueues[numQ];
  flushed_printf("Testing creation of littleQ...\n");

  for (int Log2NumCU = 0; Log2NumCU < 4; ++Log2NumCU) {
	flushed_printf("\tInit CPU_STATUS:\t");
	printBits(CPU_STATUS, numCPU);
	for (int kI = 0; kI < numQ; ++kI) {
	  bigQueues[kI] = rm->createUniqQueueLittleCore(Log2NumCU, edcl);
	  flushed_printf("\tCPU_STATUS after create littleQ(%d):\t", (int)exp2(Log2NumCU));
	  printBits(CPU_STATUS, numCPU);
	  if (bigQueues[kI] == nullptr) break;
	}

	for (int kI = 0; kI < numQ; ++kI) {
	  if (bigQueues[kI] != nullptr) {
		rm->unbindEDQueueCPU(bigQueues[kI]);
		edcl->releaseEDQueue(bigQueues[kI]);
		bigQueues[kI] = nullptr;
	  } else {
		flushed_printf("\tlittleQ(%d) ends at %d\n", (int)exp2(Log2NumCU), kI);
		break;
	  }
	}
  }
}

void testFrequencySetting() {
  flushed_printf("Testing frequency Setting...\n");
  EDQueue queues[3];
  queues[0] = rm->createUniqQueueLittleCore(2, edcl);
  queues[1] = rm->createUniqQueueBigCore(2, edcl);
  queues[2] = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  for (auto q: queues) {
	uint cur_freq = rm->getDeviceCurrentFrequency(q);
	flushed_printf("\tCurrent frequency: %d\n", cur_freq);
	for (int kI = 0;; kI++) {
	  usleep(1000000);
	  if (rm->setAssociateDeviceFrequency(q, kI) != 0) break;
	  cur_freq = rm->getDeviceCurrentFrequency(q);
	  flushed_printf("\tCurrent frequency after setting: %d\n", cur_freq);
	}
	edcl->releaseEDQueue(q);
  }
  rm->setSoCMaxFrequency();
}

void testCurrentRecording() {
  std::vector<float> current;
  int recTime = 10;
  float avgCurrent = 0;
  flushed_printf("Testing recording current (%d sec)...\n", recTime);
  if (rm->startCurrentRecording(&current) == 0) {
	sleep(recTime);
	rm->stopCurrentRecording();
  }
  flushed_printf("\tNum records: %d\n", current.size());
  for (auto &c:current) {
	avgCurrent += c;
  }
  avgCurrent /= current.size();
  flushed_printf("\tAverage current value:%f\n", avgCurrent);
}

void testFrequencySettingAndCurrentRecording() {
  flushed_printf("Testing setting frequency and recording current...\n");
  std::vector<float> current;
  int recTime = 10;
  float avgCurrent = 0;
  EDQueue queues[3];
  uint cur_freq;
  queues[0] = rm->createQueueLittleCore(2, edcl);
  queues[1] = rm->createQueueBigCore(2, edcl);
  queues[2] = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  sleep(recTime);
  for (auto q: queues) {
	for (int kI = 0;; kI++) {
	  if (rm->setAssociateDeviceFrequency(q, kI) != 0) break;
	  cur_freq = rm->getDeviceCurrentFrequency(q);
	  flushed_printf("\tCurrent frequency: %d\n", cur_freq);
	  sleep(10); // cool down 10s
	  if (rm->startCurrentRecording(&current) == 0) {
		sleep(recTime);
		rm->stopCurrentRecording();
	  }
	  flushed_printf("\t  Num records: %d. ", current.size());
	  for (auto &c:current) {
		avgCurrent += c;
	  }
	  avgCurrent /= current.size();
	  flushed_printf("Average current value:%f\n", avgCurrent);
	  current.clear();
	  avgCurrent = 0;
	}
	edcl->releaseEDQueue(q);
  }
  rm->setSoCMaxFrequency();
}

int main() {
  edcl = new EDCL();
  rm = new ResourceManager_();
  printResourceInfo();
//  coreBindingTests();
//  testRequestBigCoreQueue();
//  testRequestLittleCoreQueue();
//  testFrequencySetting();
//  testCurrentRecording();
  testFrequencySettingAndCurrentRecording();

}