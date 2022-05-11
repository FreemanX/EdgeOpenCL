#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include <iostream>

#include "CPU_utils.h"
#include "utils.h"
#include "EDCL.h"

u_int16_t CPU_STATUS = 0;
uint CPUHeteroLevel = 1; // Number of different CPU types
u_int16_t *heteroCPUMasks; // array of CPU masks for different CPU types, size=CPUHeteroLevel
uint numCPU;
void initCPUManagementEnv();
int bindEDQueueCPU(EDQueue queue);
void unbindEDQueueCPU(EDQueue queue);
void releaseCPUManagementEnv();

int main() {
  initCPUManagementEnv();
  auto edcl = new EDCL();
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
  releaseCPUManagementEnv();
  delete edcl;
  return 0;
}

void unbindEDQueueCPU(EDQueue queue) {
  flushed_printf("Unbinding queue(%d) with %d cores...\n", queue->coreStart, queue->numCU);
  flushed_printf("\tCPU_STATUS before unbind:\t");
  printBits(CPU_STATUS, numCPU);
  flushed_printf("\tqueue->numCU: %d, queue->CPUSet:\t", queue->numCU);
  printBits(queue->cpuSet, numCPU);

  for (int i = 0; i < numCPU; ++i) {
	if (IS_BIT_SET(queue->cpuSet, i) && IS_BIT_SET(CPU_STATUS, i)) {
	  flushed_printf("\tNothing to unbind, CPU_STATUS:\t");
	  printBits(CPU_STATUS, numCPU);
	  printf("Done!\n\n");
	  return;
	}
  }

  CPU_STATUS ^= queue->cpuSet;
  flushed_printf("\tCPU_STATUS:\t\t\t");
  printBits(CPU_STATUS, numCPU);
  printf("Done!\n\n");
}

int bindEDQueueCPU(EDQueue queue) {
  flushed_printf("Binding queue(%d) with %d cores...\n", queue->coreStart, queue->numCU);
  flushed_printf("\tCPU_STATUS before bind:\t\t");
  printBits(CPU_STATUS, numCPU);

  // check core availability
  for (int i = 0; i < numCPU; ++i) {
	if (IS_BIT_SET(queue->cpuSet, i) && !IS_BIT_SET(CPU_STATUS, i)) {
	  flushed_printf("\tCore binding FAIL, CPU set:\t");
	  printBits(queue->cpuSet, numCPU);
	  flushed_printf("\tCPU_STATUS:\t\t\t");
	  printBits(CPU_STATUS, numCPU);
	  printf("Done!\n\n");
	  return -1;
	}
  }

  CPU_STATUS ^= queue->cpuSet;
  flushed_printf("\tCore binding SUCCESS, CPU set:\t");
  printBits(queue->cpuSet, numCPU);
  flushed_printf("\tCPU_STATUS:\t\t\t");
  printBits(CPU_STATUS, numCPU);
  printf("Done!\n\n");
  return 0;
}

void initCPUManagementEnv() {
  flushed_printf("%s\n", __func__);
  numCPU = getNumCPUs();
  // set all CPU available, all bits to 1
  CPU_STATUS = exp2(numCPU) - 1;
  std::ifstream cpuFreqFStream;
  char cpuFreqFile[128];
  uint cpuFreqAll[numCPU]; // all cores' frequencies
  for (int i = 0; i < numCPU; ++i) {
	snprintf(cpuFreqFile, sizeof(cpuFreqFile), "/sys/devices/system/cpu/cpu%d/cpufreq/scaling_max_freq", i);
	cpuFreqFStream.open(cpuFreqFile);
	if (cpuFreqFStream.is_open()) {
	  cpuFreqFStream >> cpuFreqAll[i];
	  flushed_printf(" CPU %d freq: %d\n", i, cpuFreqAll[i]);
	  if (i > 0 && cpuFreqAll[i] != cpuFreqAll[i - 1]) CPUHeteroLevel++;
	} else {
	  flushed_printf("!Failed to Get CPU%d info from %s\n", i, cpuFreqFile);
	  exit(-1);
	}
	cpuFreqFStream.close();
  }
  heteroCPUMasks = (u_int16_t *)calloc(CPUHeteroLevel, sizeof(u_int16_t));
  int maskIdx = 0;
  for (int i = 0; i < numCPU; ++i) {
	if (i > 0 && cpuFreqAll[i] > cpuFreqAll[i - 1]) maskIdx++; // assume big core's CPU_id is also larger
	SET_BIT(heteroCPUMasks[maskIdx], i);
  }
  flushed_printf(" CPU_STATUS: ");
  printBits(CPU_STATUS, numCPU);
  flushed_printf(" CPUHeteroLevel: %d\n", CPUHeteroLevel);
  for (int i = 0; i < CPUHeteroLevel; ++i) {
	flushed_printf(" CPU Mask%d: ", i);
	printBits(heteroCPUMasks[i], numCPU);
  }
  printf("Done!\n\n");
}

void releaseCPUManagementEnv() {
  free(heteroCPUMasks);
}

