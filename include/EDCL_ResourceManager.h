#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#ifndef EDCL_INCLUDE_EDCL_RESOURCEMANAGER_H_
#define EDCL_INCLUDE_EDCL_RESOURCEMANAGER_H_
#include "EDCL.h"
#include "CPU_utils.h"
#include <dirent.h>
#include <pthread.h>
#include <unistd.h>
#include <cerrno>
#include <unordered_map>
#include "concurrentqueue.h"

// LAMBDA_X[freq][numCore]
static float LAMBDA_L[4][3] = {{0.0055034, 0.004444, 0.002963667},
															 {0.0061346, 0.004641, 0.003058667},
															 {0.0109632, 0.0045485, 0.003059},
															 {0.037222667, 0.0148975, 0.004456}};
static float LAMBDA_B[4][3] = {{0.000977667, 0.000704, 0.000442667},
															 {0.001043667, 0.000724667, 0.000453},
															 {0.004449667, 0.000708, 0.000453},
															 {0.0017764, 0.0012328, 0.000553}};
static float LAMBDA_G[4] = {0.0032248, 0.0054986, 0.0138784, 0.037222667};

static float I0[4] = {0.081474, 0.079680057, 0.075202543, 0.072300143};

struct RecordingPack {
	std::vector<float> *currentRecords;
	pthread_mutex_t *mutex_recording;
};

static int needQuit(pthread_mutex_t *mtx) {
	switch (pthread_mutex_trylock(mtx)) {
		case 0: /* if we got the lock, unlock and return 1 (true) */
			pthread_mutex_unlock(mtx);
			return 1;
		case EBUSY: /* return 0 (false) if the mutex was locked */
			return 0;
	}
	return 1;
}

static void *recordCurrent(void *arg) {
	auto pack = (RecordingPack *)arg;
	pthread_mutex_t *mx = pack->mutex_recording;
	auto currentRecords = pack->currentRecords;
	std::ifstream currentFileStream;
	float current;
	setThreadAffinity(1);
	while (!needQuit(mx)) {
		currentFileStream.open("/sys/class/power_supply/usb/input_current_now");
//	currentFileStream.open("/sys/class/power_supply/battery/current_now");
//	if (currentFileStream.is_open()) { // disable this if
		currentFileStream >> current;
		currentFileStream.close();
//	}
		currentRecords->push_back(current / 1000000.0f);
		usleep(100);
	}
	pthread_exit(nullptr);
}

using namespace moodycamel;
struct OCLKernelExecutionThread : public Trackable {
	OCLKernelExecutionThread(EDCL *edcl, EDQueue ExeQ) {
		edcl_ = edcl;
		queue = ExeQ;
		pthread_mutex_init(&runningMutex, nullptr);
		packDes = "OCLKernelExecutionThread(" + std::to_string(getTrackNum()) + ") ";
		if (ExeQ->queueType == GPUQueue) packDes += "GPU";
		else
			packDes +=
					"CPU[" + std::to_string(ExeQ->coreStart) + ", " + std::to_string(ExeQ->coreStart + ExeQ->numCU - 1) + "]Pack";
		pthread_mutex_init(&trackerMutex, nullptr);
		pthread_mutex_init(&haltMutex, nullptr);
	}

	~OCLKernelExecutionThread() {
		pthread_mutex_destroy(&runningMutex);
		pthread_mutex_destroy(&haltMutex);
		pthread_mutex_destroy(&trackerMutex);
	}

	void updateWorkTracker(int kernelID) {
		pthread_mutex_lock(&trackerMutex);
		if (workTracker.count(kernelID) > 0) {
			workTracker[kernelID] += 1;
		} else {
			workTracker[kernelID] = 1;
		}
//	debug_printf(__func__, "[kernel %d] executed: %d, total: %d.\n",
//				 kernelID, workTracker[kernelID], edcl_->getNumKernelSlices(kernelID));
		pthread_mutex_unlock(&trackerMutex);
	}

	uint getWorkTrackerRecord(int kernelID) {
		u_int numSlice = 0;
		if (pthread_mutex_trylock(&trackerMutex) == 0) {
			if (workTracker.count(kernelID) > 0) numSlice = workTracker[kernelID];
//	debug_printf(__func__, "\n");
			pthread_mutex_unlock(&trackerMutex);
		}
		return numSlice;
	}

	bool isLaunched = false;
	std::string packDes;
	pthread_t pthread{};
	pthread_mutex_t runningMutex{};
	pthread_mutex_t haltMutex{};
	EDQueue queue;
	EDCL *edcl_;
	ConcurrentQueue<EDKernel> kernels;
	pthread_mutex_t trackerMutex{};
	std::unordered_map<int, uint> workTracker; // count number of slices of each kernel have been executed
	double activeTime = 0;
	std::vector<EDKernel> executedKernels;

	void printActiveTime() const {
		flushed_printf("%s: active time %lf\n", packDes.c_str(), activeTime);
	}
};


// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()

static void *OCLKernelExecutingThreadFuc(void *arg) {
	auto *pack = (OCLKernelExecutionThread *)arg;
	pthread_mutex_t *mx = &pack->runningMutex;
	EDQueue queue = pack->queue;
	auto edcl = pack->edcl_;
	auto kernels = &pack->kernels;
	auto executedKernels = &pack->executedKernels;
	pack->activeTime = getCurrentTime();
//  int kernelExecuted = 0;
//  debug_printf(__func__, "%s thread started! num kernels: %d\n", pack->packDes.c_str(), kernels->size_approx());
	while (kernels->size_approx() > 0 || !needQuit(mx)) {
		if (needQuit(&pack->haltMutex)) break;
		EDKernel kernel;
		if (kernels->try_dequeue(kernel)) {
//	  ++kernelExecuted;
//	  kernel->printInfo();
//	  debug_printf(__func__, "%s get 1 kernel!, total:%d\n", pack->packDes.c_str(), ++kernelExecuted);
			try {
				edcl->executeKernel(queue,
														kernel,
														kernel->num_events_in_wait_list,
														kernel->event_wait_list,
														&kernel->kernelEvent);
			} catch (std::runtime_error &e) {
				err_printf(__func__, __LINE__, "Pack %s: %s\n", pack->packDes.c_str(), e.what());
				err_printf(__func__, __LINE__, "Kernel %s info: offset(%d, %d, %d), global(%d, %d, %d), local(%d, %d, %d) \n",
									 kernel->name.c_str(),
									 CHOOSE(kernel->offset != nullptr, kernel->offset[1], 0),
									 CHOOSE(kernel->offset != nullptr && kernel->workDim > 1, kernel->offset[1], 0),
									 CHOOSE(kernel->offset != nullptr && kernel->workDim > 2, kernel->offset[2], 0),
									 kernel->globalSize[0],
									 CHOOSE(kernel->workDim > 1, kernel->globalSize[1], 0),
									 CHOOSE(kernel->workDim > 2, kernel->globalSize[2], 0),
									 kernel->localSize[0],
									 CHOOSE(kernel->workDim > 1, kernel->localSize[1], 0),
									 CHOOSE(kernel->workDim > 2, kernel->localSize[2], 0)
				);
			}
			pack->updateWorkTracker(kernel->kernelID);
			executedKernels->push_back(kernel);
		}
	}
//  debug_printf(__func__, "%s thread end! num kernels left: %d\n", pack->packDes.c_str(), kernels->size_approx());
//  debug_printf(__func__, "%s break out of while loop, executed kernels: %d\n", pack->packDes.c_str(), kernelExecuted);
	for (auto &k: *executedKernels) {
//	debug_printf(__func__, "%s waiting for kernel %s...\n",
//				 pack->packDes.c_str(), k->name.c_str());
		k->waitForKernelEvent();
//	debug_printf(__func__, "%s kernel %s execution time: %lf\n",
//				 pack->packDes.c_str(), k->name.c_str(), k->getKernelEventTime());
	}
	pack->activeTime = getCurrentTime() - pack->activeTime;
	pthread_exit(nullptr);
}

typedef class ResourceManager_ *ResourceManager;
// TODO: Make all functions thread safe
class ResourceManager_ {
 public:
	ResourceManager_() {
		initCPUManagementEnv();
	}

	~ResourceManager_() {
		setSoCMaxFrequency();
		free(heteroCPUMasks);
	}

	std::unordered_map<int, EDQueue> uniqueLQueues;
	std::unordered_map<int, EDQueue> uniqueBQueues;
	int bindEDQueueCPU(EDQueue queue) {
		// bind check if the core has been used
		for (int i = 0; i < numCPU; ++i)
			if (IS_BIT_SET(queue->cpuSet, i) && !IS_BIT_SET(AVAILABLE_CPU, i)) return -1;
		AVAILABLE_CPU ^= queue->cpuSet;
		if (queue->coreStart < numLittleCore) uniqueLQueues[queue->getTrackNum()] = queue;
		else uniqueBQueues[queue->getTrackNum()] = queue;
		return 0;
	}

	void releaseNoBindUniqQueue(EDQueue queue, EDCL *edcl) {
		if (queue == nullptr) return;
		if (uniqueLQueues.count(queue->getTrackNum()) > 0) uniqueLQueues.erase(queue->getTrackNum());
		else if (uniqueBQueues.count(queue->getTrackNum())) uniqueBQueues.erase(queue->getTrackNum());
		edcl->releaseEDQueue(queue);
	}

	void unbindEDQueueCPU(EDQueue queue) {
		// unbind check if the cores aren't bind before, then :
		for (int i = 0; i < numCPU; ++i)
			if (IS_BIT_SET(queue->cpuSet, i) && IS_BIT_SET(AVAILABLE_CPU, i)) return;
		if (queue->coreStart < numLittleCore) uniqueLQueues.erase(queue->getTrackNum());
		else uniqueBQueues.erase(queue->getTrackNum());
		AVAILABLE_CPU ^= queue->cpuSet;
	}

	uint getCPUHeteroLevel() const { return this->numHeteroCores.size(); }
	uint getNumCPU() const { return this->numCPU; }
	std::vector<uint> getHeteroCoreNums() { return this->numHeteroCores; }
	u_int16_t getAVAILABLE_CPU() const { return AVAILABLE_CPU; }
	uint getNumBigCores() const { return numBigCore; }
	uint getNumLittleCores() const { return numLittleCore; }

	uint getCPUCoreHeteroLevel(int cpuID) {
		uint coreCnt = 0;
		for (int i = 0; i < getCPUHeteroLevel(); ++i) {
			coreCnt += numHeteroCores[i];
			if (cpuID < coreCnt) return i;
		}
		return 65536;

	}

	uint getEDQueueHeteroLevel(EDQueue queue) {
		if (queue->queueType == GPUQueue) return 0;
		if (queue->numCU == numCPU) return 65535;
		return getCPUCoreHeteroLevel(queue->coreStart);
	}

	EDQueue createQueueLittleCore(uint Log2NumCU, EDCL *edcl) const {
		uint numCore = exp2(Log2NumCU);
		if (numCore > numLittleCore) return nullptr;
		int fcid = (int)(numLittleCore - numCore);
		uint SDIndex = fcid / numCore;
		EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
		return queue;
	}

	EDQueue createQueueBigCore(uint Log2NumCU, EDCL *edcl) const {
		uint numCore = exp2(Log2NumCU);
		if (numCore > numBigCore) return nullptr;
		uint fcid = numCPU - numCore;
		uint SDIndex = fcid / numCore;
		EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
		return queue;
	}

	EDQueue createUniqQueueLittleCore(uint Log2NumCU, EDCL *edcl) {
		uint numCore = exp2(Log2NumCU);
		if (numCore > numLittleCore) return nullptr;
		for (int fcid = (int)(numLittleCore - numCore); fcid >= 0; fcid -= numCore) {
			uint cpu_available = 0;
			for (int cpuID = fcid; cpuID < fcid + numCore; ++cpuID) {
				cpu_available += IS_BIT_SET(AVAILABLE_CPU, cpuID); // not available
			}
			if (cpu_available != numCore) continue;

			uint SDIndex = fcid / numCore;
//	  flushed_printf("fcid: %d, numCore: %d, SDIndex: %d\n", fcid, numCore, SDIndex);
			EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
			bindEDQueueCPU(queue);
			return queue;
		}
		return nullptr;
	}

	EDQueue createUniqQueueBigCore(uint Log2NumCU, EDCL *edcl) {
		/*
		 * Number of Big cores(NBC): numBigCore
		 * Possible Log2NumCU(PL) = 0, 1, 2.
		 * Number of cores applied(NC): exp2(PL). NC < NBC
		 *
		 * CPU ids: 4, 5, 6, 7.
		 * First CPU ID(FCID): NumCPU(8) - NBC(4) = 4, NumCPU(8) - NBC(2) = 6
		 * First SDIndex: FCID/numCores
		 * */
		uint numCore = exp2(Log2NumCU);
		if (numCore > numBigCore) return nullptr;
//	for (uint fcid = numCPU - numBigCore; fcid < numCPU; fcid += numCore) {
		for (uint fcid = numCPU - numCore; fcid >= numCPU - numBigCore; fcid -= numCore) {
			uint cpu_available = 0;
			for (int cpuID = fcid; cpuID < fcid + numCore; ++cpuID) {
				cpu_available += IS_BIT_SET(AVAILABLE_CPU, cpuID); // not available
			}
			if (cpu_available != numCore) continue;
			uint SDIndex = fcid / numCore;
			EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
			bindEDQueueCPU(queue);
			return queue;
		}
		return nullptr;
	}

	int getUniqQueueUsedLittleCores() {
		int numCores = 0;
		for (auto &qp : uniqueLQueues) {
			auto q = qp.second;
			if (q->queueType == SDQueue) numCores += q->numCU;
			else uniqueLQueues.erase(qp.first);
//			numCores += q->numCU;
		}
		if (numCores > 4) {
			err_printf(__func__, __LINE__, "numCores = %d!!, uniqueLQueues size: %d\n",
								 numCores, uniqueLQueues.size());
			for (auto &qp : uniqueLQueues) {
				auto q = qp.second;
				err_printf(__func__, __LINE__, "LQueue(%d) type: %d, coreStart: %d, numCU: %d\n",
									 q->getTrackNum(), q->queueType, q->coreStart, q->numCU);
			}
		}
		assert(numCores <= 4);
		return numCores;
	}

	int getUniqQueueUsedBigCores() {
		int numCores = 0;
		for (auto &qp : uniqueBQueues) {
			auto q = qp.second;
			if (q->queueType == SDQueue) numCores += q->numCU;
			else uniqueBQueues.erase(qp.first);
		}
		if (numCores > 4) {
			err_printf(__func__, __LINE__, "numCores = %d!!, uniqueBQueues size: %d\n",
								 numCores, uniqueBQueues.size());
			for (auto &qp : uniqueBQueues) {
				auto q = qp.second;
				err_printf(__func__, __LINE__, "BQueue(%d) type: %d, coreStart: %d, numCU: %d\n",
									 q->getTrackNum(), q->queueType, q->coreStart, q->numCU);
			}
		}
		assert(numCores <= 4);
		return numCores;
	}

	EDQueue createUniqQueueLittleCore_noBind(uint Log2NumCU, EDCL *edcl) {
		int numCore = exp2(Log2NumCU);
		if (numCore > numLittleCore - getUniqQueueUsedLittleCores()) return nullptr;
		uint SDIndex = (numLittleCore - numCore) / numCore;
		EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
		uniqueLQueues[queue->getTrackNum()] = queue;
		return queue;
	}

	EDQueue createUniqQueueBigCore_noBind(uint Log2NumCU, EDCL *edcl) {
		int numCore = exp2(Log2NumCU);
		if (numCore > numBigCore - getUniqQueueUsedBigCores()) return nullptr;
		uint SDIndex = (numCPU - numCore) / numCore;
		EDQueue queue = edcl->createSubDeviceExeQueue(Log2NumCU, SDIndex);
		uniqueBQueues[queue->getTrackNum()] = queue;
		return queue;
	}

	void releaseUniqQueue(EDQueue Q, EDCL *edcl) {
		if (Q == nullptr) return;
//	Q->finish();
		if (Q->queueType == SDQueue) unbindEDQueueCPU(Q);
		edcl->releaseEDQueue(Q);
	}

	size_t getNumAvailableFrequencies(EDQueue queue) {
		if (queue->queueType == GPUQueue) return GPUFrequencies.size();
		if (queue->numCU == numCPU) return 65536;
		return heteroCPUFrequencies[getEDQueueHeteroLevel(queue)].size();
	}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "readability-convert-member-functions-to-static"
	uint getGPUCurrentFrequency() const {
		std::ifstream freqFStream;
		uint freq;
		freqFStream.open("/sys/class/kgsl/kgsl-3d0/gpuclk");
		if (freqFStream.is_open()) {
			freqFStream >> freq;
			freqFStream.close();
			return freq / 1000000;
		} else {
			flushed_printf("Opening GPU frequency file failed!\n");
			return 65536;
		}
	}
#pragma clang diagnostic pop

	uint getCPUCurrentFrequency(int cpuID) const {
		if (cpuID >= numCPU) {
			err_printf(__func__, __LINE__, "CPU core %d doesn't exist!\n", cpuID);
			return 65535;
		}
		std::ifstream freqFStream;
		uint freq;
		char cpuFreqFile[128];
		snprintf(cpuFreqFile,
						 sizeof(cpuFreqFile),
						 "/sys/devices/system/cpu/cpu%d/cpufreq/scaling_cur_freq",
						 cpuID);
		freqFStream.open(cpuFreqFile);
		if (freqFStream.is_open()) {
			freqFStream >> freq;
			freqFStream.close();
			return freq;
		} else {
			err_printf(__func__, __LINE__, "Opening CPU frequency file failed!\n");
			return 65536;
		}
	}

	uint getDeviceCurrentFrequency(EDQueue queue) const {
		if (queue->queueType == GPUQueue) return getGPUCurrentFrequency();
		if (queue->numCU == numCPU) return 65536; // TODO: this should check if CUs contain hetero cores
		return getCPUCurrentFrequency(queue->coreStart);
	}

	int setGPUFrequency(int freqIdx) {
		if (freqIdx >= GPUFrequencies.size()) return -1;
		char cmd_setMax[128];
		char cmd_setMin[128];
		uint curFreq = getGPUCurrentFrequency();
		uint targetFreq = GPUFrequencies[freqIdx];
		snprintf(cmd_setMax, sizeof(cmd_setMax),
						 "echo %d > /sys/class/kgsl/kgsl-3d0/max_clock_mhz", targetFreq);
		snprintf(cmd_setMin, sizeof(cmd_setMin),
						 "echo %d > /sys/class/kgsl/kgsl-3d0/min_clock_mhz", targetFreq);
		return setFrequency(targetFreq, curFreq, cmd_setMin, cmd_setMax);
	}

	int setCPUFrequency(uint cpuID, uint freqIdx) {
		if (cpuID >= numCPU) return -1;
		uint h_level = getCPUCoreHeteroLevel(cpuID);
		if (freqIdx >= heteroCPUFrequencies[h_level].size()) return -2;
		char cmd_setMax[128];
		char cmd_setMin[128];
		uint curFreq = getCPUCurrentFrequency(cpuID);
		uint targetFreq = heteroCPUFrequencies[h_level].at(freqIdx);
		snprintf(cmd_setMax,
						 sizeof(cmd_setMax),
						 "echo %d > /sys/devices/system/cpu/cpu%d/cpufreq/scaling_max_freq",
						 targetFreq, cpuID);
		snprintf(cmd_setMin,
						 sizeof(cmd_setMin),
						 "echo %d > /sys/devices/system/cpu/cpu%d/cpufreq/scaling_min_freq",
						 targetFreq, cpuID);
		return setFrequency(targetFreq, curFreq, cmd_setMin, cmd_setMax);
	}

	void setSoCMaxFrequency() {
		for (int i = 0; i < numCPU; ++i) setCPUFrequency(i, 0);
		setGPUFrequency(0);
	}

	int setAssociateDeviceFrequency(EDQueue queue, int freqIdx) {
		if (queue->queueType == GPUQueue) return setGPUFrequency(freqIdx);
		if (queue->numCU == numCPU) return -2;
		return setCPUFrequency(queue->coreStart, freqIdx);
	}

	std::vector<std::vector<uint>> getHeteroCPUFrequencies() { return heteroCPUFrequencies; }
	std::vector<uint> getGPUFrequencies() { return GPUFrequencies; }

	int startCurrentRecording(std::vector<float> &current) {
		if (pthread_mutex_trylock(&mutex_createThread_) == 0) {
			pthread_mutex_lock(&mutex_recording_);
			auto *rp = new RecordingPack;
			rp->currentRecords = &current;
			rp->mutex_recording = &mutex_recording_;
			pthread_create(&currentRecordingThread_, nullptr, recordCurrent, rp);
			return 0;
		} else {
			flushed_printf("Current recording has been launched!\n");
			return -1;
		}
	}

	int startCurrentRecording(std::vector<float> *current) {
		if (pthread_mutex_trylock(&mutex_createThread_) == 0) {
			pthread_mutex_lock(&mutex_recording_);
			auto *rp = new RecordingPack;
			rp->currentRecords = current;
			rp->mutex_recording = &mutex_recording_;
			pthread_create(&currentRecordingThread_, nullptr, recordCurrent, rp);
			return 0;
		} else {
			flushed_printf("Current recording has been launched!\n");
			return -1;
		}
	}

	void stopCurrentRecording() {
		pthread_mutex_unlock(&mutex_recording_);
		pthread_join(currentRecordingThread_, nullptr);
		pthread_mutex_unlock(&mutex_createThread_);
	}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "readability-convert-member-functions-to-static"
	void startOCLKernelExecutionThread(OCLKernelExecutionThread *exeThread) {
		if (exeThread->isLaunched) return;
		OCLKernelExecutionThread *thread = exeThread;
		pthread_mutex_lock(&thread->runningMutex);
		pthread_mutex_lock(&thread->haltMutex);
		int rc;
		if ((rc = pthread_create(&thread->pthread, nullptr, OCLKernelExecutingThreadFuc, thread))) {
			fprintf(stderr, "[%s]error: pthread_create, rc: %d\n", __func__, rc);
			exit(EXIT_FAILURE);
		} else exeThread->isLaunched = true;
	}

	void haltRunningThread(OCLKernelExecutionThread *thread) {
		pthread_mutex_unlock(&thread->haltMutex);
		pthread_mutex_unlock(&thread->runningMutex);
	}

	void unlockRunningThread(OCLKernelExecutionThread *thread) {
		pthread_mutex_unlock(&thread->runningMutex);
	}

#pragma ide diagnostic ignored "bugprone-unused-return-value"
	void joinOCLKernelExecutionThread(OCLKernelExecutionThread *executionThread) {
		if (!executionThread->isLaunched) return;
		pthread_join(executionThread->pthread, nullptr);
		executionThread->isLaunched = false;
		pthread_mutex_trylock(&executionThread->haltMutex); // keep it consistent with runningMutex
		pthread_mutex_unlock(&executionThread->haltMutex);
	}
#pragma clang diagnostic pop

 private:
	u_int16_t AVAILABLE_CPU = 0; // mask
	std::vector<uint> numHeteroCores; // stores number of each type of cores, size of vector is CPUHeteroLevel
	std::vector<std::vector<uint>> heteroCPUFrequencies;
	std::vector<uint> GPUFrequencies;
	uint numBigCore = 0;
	uint numLittleCore = 0;
	u_int16_t *heteroCPUMasks{}; // array of CPU masks for different CPU types, size=CPUHeteroLevel
	uint numCPU{};
	pthread_t currentRecordingThread_{};
	pthread_mutex_t mutex_recording_{}, mutex_createThread_{};

	void initCPUManagementEnv() {
		pthread_mutex_init(&mutex_createThread_, nullptr);
		pthread_mutex_init(&mutex_recording_, nullptr);
		numCPU = getNumCPUs();
		// set all CPU available, all bits to 1
		AVAILABLE_CPU = exp2(numCPU) - 1;
		// Distinguish Big, Little core via CPU Frequency
		std::ifstream cpuFreqFStream;
		char cpuFreqFile[128];
		uint cpuFreqAll[numCPU]; // all cores' max frequencies
		uint numCores = 0;
		for (int i = 0; i < numCPU; ++i) {
			snprintf(cpuFreqFile, sizeof(cpuFreqFile), "/sys/devices/system/cpu/cpu%d/cpufreq/cpuinfo_max_freq", i);
			cpuFreqFStream.open(cpuFreqFile);
			if (cpuFreqFStream.is_open()) {
				cpuFreqFStream >> cpuFreqAll[i];
				if (i > 0 && cpuFreqAll[i] != cpuFreqAll[i - 1]) { // if encounters different max freq -> new core type
					numHeteroCores.push_back(numCores);
					numCores = 1;
				} else numCores++;

			} else {
				flushed_printf("!Failed to Get CPU%d freq from %s\n", i, cpuFreqFile);
				exit(-1);
			}
			cpuFreqFStream.close();
		}
		numHeteroCores.push_back(numCores);
		heteroCPUMasks = (u_int16_t *)calloc(numHeteroCores.size(), sizeof(u_int16_t));
		int maskIdx = 0;
		for (int i = 0; i < numCPU; ++i) {
			if (i > 0 && cpuFreqAll[i] > cpuFreqAll[i - 1]) maskIdx++; // assume big core's CPU_id is also larger
			SET_BIT(heteroCPUMasks[maskIdx], i);
		}
		numLittleCore = numHeteroCores[0];
		for (int i = 1; i < getCPUHeteroLevel(); ++i) {
			numBigCore += numHeteroCores[i];
		}
		loadGPUFrequencies();
		loadHeteroCPUFrequencies();
	}

//  std::vector<std::vector<uint>> heteroCPUFrequencies;
	void loadHeteroCPUFrequencies() {
		// path: /sys/devices/system/cpu/cpufreq/policy*/scaling_available_frequencies
		// path: /sys/devices/system/cpu/cpufreq/policy4/scaling_available_frequencies
		// max freq path: /sys/devices/system/cpu/cpufreq/policy4/cpuinfo_max_freq
		std::string cpufreqPath = "/sys/devices/system/cpu/cpufreq/";
		std::string cpuAvailableFreqFileName = "/scaling_available_frequencies";
		std::string cpuMaxFreqFileName = "/cpuinfo_max_freq";
		struct dirent *de;  // Pointer for directory entry
		// opendir() returns a pointer of DIR type.
		DIR *dr = opendir(cpufreqPath.data());

		if (dr == nullptr) {
			printf("Could not open /sys/devices/system/cpu/cpufreq/ directory");
			return;
		}

		// Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html
		std::vector<uint> *freqData;
		while ((de = readdir(dr)) != nullptr) {
			if (std::string(de->d_name).find('.') == std::string::npos) {
				std::string policyDir = cpufreqPath + de->d_name;
				freqData = load_uint_from_file(policyDir + cpuAvailableFreqFileName);
				if (freqData != nullptr) {
					std::reverse(freqData->begin(), freqData->end());
					std::vector<uint> CPUFreqs;
					for (int kI = 1; kI < freqData->size() && kI <= 10; kI += 4) {
						uint freq = freqData->at(kI);
						CPUFreqs.push_back(freq);
					}
					delete freqData;
					freqData = load_uint_from_file(policyDir + cpuMaxFreqFileName);
					CPUFreqs.push_back(freqData->at(0));
					reverseSortUniqueVector(CPUFreqs);
					heteroCPUFrequencies.push_back(CPUFreqs);
					delete freqData;
				} else {
					err_printf(__func__, __LINE__, "Can't find available CPU frequencies file from \"%s\"!\n", policyDir.data());
				}
			}
		}
		// little core has lower hetero level num. readdir get dir from higher policy num
		std::reverse(heteroCPUFrequencies.begin(), heteroCPUFrequencies.end());
		closedir(dr);
	}

	static void reverseSortUniqueVector(std::vector<uint> &vec) {
		// remove duplicate
		vec.erase(unique(vec.begin(), vec.end()), vec.end());
		// sort
		std::sort(vec.begin(), vec.end());
		// reverse order
		std::reverse(vec.begin(), vec.end());
	}

//  std::vector<uint> GPUFrequencies;
	void loadGPUFrequencies() {
		// path: /sys/class/kgsl/kgsl-3d0/gpu_available_frequencies
		std::vector<uint> *data = load_uint_from_file("/sys/class/kgsl/kgsl-3d0/gpu_available_frequencies");
		if (data != nullptr) {
			for (int kI = 0; kI < data->size(); kI += 2) {
				uint freq = data->at(kI) / 1000000;
				GPUFrequencies.push_back(freq);
			}
			reverseSortUniqueVector(GPUFrequencies);
		} else {
			flushed_printf("Can't find available GPU frequencies file!\n");
		}
		delete data;
	}

	static std::vector<uint> *load_uint_from_file(std::string &&filePath) {
		std::ifstream infile(filePath);
		std::vector<uint> *vec = nullptr;
		if (infile) {
			vec = new std::vector<uint>((std::istream_iterator<uint>(infile)), std::istream_iterator<uint>());
			infile.close();
		}
		return vec;
	}

	static int setFrequency(uint targetFreq, uint curFreq, char *cmd_setMin, char *cmd_setMax) {
		int status = 0;
		if (targetFreq < curFreq) { // decrease frequency
			status += std::system(cmd_setMin);
			status += std::system(cmd_setMax);
		} else if (targetFreq > curFreq) { // increase frequency
			status += std::system(cmd_setMax);
			status += std::system(cmd_setMin);
		}
		return status;
	}

};

#undef globalSize
#undef localSize
#undef offset
#endif //EDCL_INCLUDE_EDCL_RESOURCEMANAGER_H_