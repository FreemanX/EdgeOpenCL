#pragma ide diagnostic ignored "OCDFAInspection"
#ifndef EDCL_INCLUDE_STG_CA2020_H_
#define EDCL_INCLUDE_STG_CA2020_H_
#include "EDCL_Strategy.h"
#include <unordered_map>
#include "csv.h"
#include <memory>
#include <queue>
#include <utility>

static int BENCHMARK_KERNEL_SLICE_LOOP_NUM = 4; // number of executions for generating a DP
static double SET_DP_TIME = 2.0; // time in sec a kernel should run when generating a DP
#define TDP_FILE_NAME "PolyTDP_New.csv"
#define TDP_FILE_DIS_NAME "PolyTDP_New_distill.csv"

enum Preference { PERFORMANCE, ENERGY, PRODUCT };

// Design Points
struct DP {
	std::string DP_ID;
	float Prf{}; // performance
	// Energy consumptions
//  float EC_L{};
//  float EC_B{};
//  float EC_G{};
// EC can be calculated by frequencies and current(total).
	float EC{};
//  std::vector<float> currentLog;
	int NC{}; // number of chunks. Stored in form of power2 factor
	std::vector<float> CK_ratio{0, 0, 0}; // portion[0,1] of chunks to be executed on each devices.
	std::vector<size_t> CK{0, 0, 0}; // number of chunks to be executed on each devices.
	int nb = -1; // -1 indicates not using this core. Log2Num: 0,1,2
	int nl = -1; // -1 indicates not using this core. Log2Num: 0,1,2
	int ng = -1; // -1 indicates not using this core. Log2Num: 0,1. Indicating using GPU or not
	// Frequencies -> related to energy consumption
	int fb = 0;
	int fl = 0;
	int fg = 0;
	bool operator<(const DP &dp) const {
		return (Prf < dp.Prf);
	}
	bool operator>(const DP &dp) const {
		return (Prf > dp.Prf);
	}
	// runtime
	bool bindWithKernel = false;
	void printDP() {
		flushed_printf("%s,", DP_ID.c_str());
		flushed_printf("%f,", Prf);
		flushed_printf("%f,", EC);
		flushed_printf("%d,", NC);

		flushed_printf("[ ");
		for (auto ratio:CK_ratio) flushed_printf("%f ", ratio);
		flushed_printf("],");

		flushed_printf("[ ");
		for (auto ck:CK) flushed_printf("%d ", ck);
		flushed_printf("],");

		flushed_printf("%d,", nl);
		flushed_printf("%d,", fl);
		flushed_printf("%d,", nb);
		flushed_printf("%d,", fb);
		flushed_printf("%d,", ng);
		flushed_printf("%d\n", fg);
	}
};

using namespace moodycamel;
struct KernelExecutionThread : public Trackable {
	KernelExecutionThread(EDCL *edcl, ResourceManager RM, DP &dp, EDKernel Kernel) {
		edcl_ = edcl;
		this->dp = dp;
		kernel = Kernel;
		rm = RM;
		queues.resize(3);
		queues[0] = nullptr;
		queues[1] = nullptr;
		queues[2] = nullptr;
		threadDescription = "KernelExecutionThread(" + std::to_string(getTrackNum()) + ") " + Kernel->name + "["
				+ std::to_string(Kernel->kernelID) + "]";
		pthread_mutex_init(&trackerMutex, nullptr);
		pthread_mutex_init(&runningMutex, nullptr);
	}

	void setDP(DP &dp_new) { this->dp = dp_new; }

	bool prepareExecution() {
		queues[0] = CHOOSE(dp.nl >= 0, rm->createUniqQueueLittleCore_noBind(dp.nl, edcl_), nullptr);
		if (dp.nl >= 0 && queues[0] == nullptr) {
			debug_printf(__func__, "%s requesting %d Little cores failed!\n", threadDescription.c_str(), pow2(dp.nl));
			return false;
		}
		queues[1] = CHOOSE(dp.nb >= 0, rm->createUniqQueueBigCore_noBind(dp.nb, edcl_), nullptr);
		if (dp.nb >= 0 && queues[1] == nullptr) {
			releaseQueues();
			debug_printf(__func__, "%s requesting %d Big cores failed!\n", threadDescription.c_str(), pow2(dp.nb));
			return false;
		}
		queues[2] = CHOOSE(dp.ng >= 0, edcl_->createDeviceCmdQueueProfilingEnabled(GPUQueue), nullptr);
		if (chunks == nullptr) chunks = edcl_->slicingKernelPow2_extra_equal(kernel, dp.NC);
		if (chunks == nullptr) {
			releaseQueues();
			debug_printf(__func__, "%s can't be sliced with slicing factor: %d\n", threadDescription.c_str(), dp.NC);
			return false;
		}
		return true;
	}

	void releaseQueues() {
		if (queues[0] != nullptr) {
			rm->releaseNoBindUniqQueue(queues[0], edcl_);
			queues[0] = nullptr;
		}
		if (queues[1] != nullptr) {
			rm->releaseNoBindUniqQueue(queues[1], edcl_);
			queues[1] = nullptr;
		}
		if (queues[2] != nullptr) {
			edcl_->releaseEDQueue(queues[2]);
			queues[2] = nullptr;
		}
	}

	~KernelExecutionThread() {
		if (queues[0] != nullptr) rm->releaseNoBindUniqQueue(queues[0], edcl_);
		if (queues[1] != nullptr) rm->releaseNoBindUniqQueue(queues[1], edcl_);
		if (queues[2] != nullptr) edcl_->releaseEDQueue(queues[2]);
		pthread_mutex_destroy(&runningMutex);
		pthread_mutex_destroy(&trackerMutex);
		if (chunks != nullptr) {
			// TODO: update kernelSliceTracker
			for (auto &k : chunks->kernels) {
				delete k;
			}
		}
	}

	void addExecutedKernel(EDKernel k) {
		pthread_mutex_lock(&trackerMutex);
		executedChunks.push_back(k);
		pthread_mutex_unlock(&trackerMutex);
	}

	uint getNumExecutedKernels() {
		u_int numSlice = 0;
		if (pthread_mutex_trylock(&trackerMutex) == 0) {
			numSlice = executedChunks.size();
			pthread_mutex_unlock(&trackerMutex);
		}
		return numSlice;
	}

	bool allSlicesEnqueued() {
		return getNumExecutedKernels() == chunks->kernels.size();
	}

	bool kernelCompleted() {
		return allSlicesEnqueued() && !isLaunched;
	}

	ResourceManager rm;
	bool isLaunched = false;
	std::string threadDescription;
	pthread_t pthread{};
	pthread_mutex_t runningMutex{};
	std::vector<EDQueue> queues;
	AtomicKernelSet chunks = nullptr; // to be executed
	EDCL *edcl_;
//	DPKernelPair *dp_kernel_pair;
	DP dp;
	EDKernel kernel;
	pthread_mutex_t trackerMutex{};
	std::vector<EDKernel> executedChunks; // record slices have been executed
	double activeTime = 0;
	double CT = 0;

	void initThread() {
		executedChunks.clear();
		activeTime = 0;
		CT = 0;
		pthread_mutex_trylock(&runningMutex); // NOLINT(bugprone-unused-return-value)
		pthread_mutex_unlock(&runningMutex);
	}

	void printActiveTime() const {
		printf("%s: active time %lf, Cumulate chunk time: %f\n",
					 threadDescription.c_str(), activeTime, chunks->atomicExeTime);
	}
};

#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()
#define ONE_CHUNK_AT_A_TIME
static void *KernelExecutionThreadFunction(void *arg) {
	auto thread = (KernelExecutionThread *)arg;
//	debug_printf(__func__, "%s starts!\n", thread->threadDescription.c_str());
	auto kernel = thread->kernel;
	auto edcl = thread->edcl_;
	auto rm = thread->rm;
	auto chunks = thread->chunks;

	double activeTime = getCurrentTime();

	auto dp = &thread->dp;
	if (dp->nl >= 0) rm->setCPUFrequency(0, dp->fl);
	if (dp->nb >= 0) rm->setCPUFrequency(rm->getNumCPU() - 1, dp->fb);
	if (dp->ng >= 0) rm->setGPUFrequency(dp->fg);
	std::vector<EDKernel> runningKernels;
	// Distribute kernels
	int largestDim = 0;
	int lastUsedDevice = 0;
	size_t totalNumWG;
	{
		auto k = chunks->kernels[0];
		for (int kJ = 0; kJ < k->workDim; ++kJ) {
			if (k->globalSize[kJ] / k->localSize[kJ] > k->globalSize[largestDim] / k->localSize[largestDim])
				largestDim = kJ;
		}
		totalNumWG = k->globalSize[largestDim] / k->localSize[largestDim];
		for (int kI = 0; kI < 3; ++kI) {
			if (dp->CK_ratio[kI] > 0) lastUsedDevice = kI;
		}
	}

//	debug_printf(__func__, "largestDim: %d, lastUsedDevice: %d, totalNumWG: %d\n",
//							 largestDim, lastUsedDevice, totalNumWG);

	std::unique_ptr<size_t[]> k_offset = std::make_unique<size_t[]>(kernel->workDim);
	std::unique_ptr<size_t[]> k_global = std::make_unique<size_t[]>(kernel->workDim);
	double singleChunkT;
	// resume kernel execution if it's been started before
	for (int i = thread->executedChunks.size(); i < chunks->kernels.size(); i++) {
		auto k = chunks->kernels[i];
		if (needQuit(&thread->runningMutex)) break;
		singleChunkT = getCurrentTime();
		memset(k_offset.get(), 0, k->workDim * sizeof(size_t));
		memset(k_global.get(), 0, k->workDim * sizeof(size_t));
//		debug_printf(__func__, "Handling slice[%d], total: %d\n", i, chunks->kernels.size());
		for (int kI = 0; kI < 3; ++kI) { // kI: loop through processors
			if (dp->CK_ratio[kI] > 0 && thread->queues[kI] != nullptr) {
				EDKernel processorSlice = nullptr;
				for (int kJ = 0; kJ < k->workDim; ++kJ) { //kJ: loop through workDim
					if (kJ == largestDim) {
						size_t numSliceWG = floor((float)totalNumWG * dp->CK_ratio[kI]);
						if (numSliceWG > 0 || (kI == lastUsedDevice && k_offset[kJ] < k->globalSize[kJ])) {
							processorSlice = edcl->createKernelCopy(k);
							k_global[kJ] = numSliceWG * k->localSize[kJ];
						}
					} else k_global[kJ] = k->globalSize[kJ];
				}
				if (thread->queues[kI] != nullptr && processorSlice == nullptr) {
					if (kI == 2) {
						edcl->releaseEDQueue(thread->queues[kI]);
					} else {
						rm->releaseNoBindUniqQueue(thread->queues[kI], edcl);
					}
					thread->queues[kI] = nullptr;
				}
				if (processorSlice == nullptr) continue;
				for (int kJ = 0; kJ < k->workDim; ++kJ) { // loop through work-dim, check if all WIs are included
					processorSlice->offset[kJ] = k_offset[kJ] + k->getOffset()[kJ];
					if (kI == lastUsedDevice && k->globalSize[kJ] - k_offset[kJ] > 0)
						processorSlice->globalSize[kJ] = k->globalSize[kJ] - k_offset[kJ];
					else
						processorSlice->globalSize[kJ] = k_global[kJ];
					if (kJ == largestDim) {
						k_offset[kJ] += k_global[kJ];
					}
				}
				dp->CK[kI] = processorSlice->globalSize[largestDim] / processorSlice->localSize[largestDim];
				assert(dp->CK[kI] > 0);
//				processorSlice->printInfo();
				edcl->confirmExeEnv(thread->queues[kI], processorSlice);
				try {
					edcl->executeKernel(thread->queues[kI], processorSlice, 0, nullptr, &processorSlice->kernelEvent);
					runningKernels.push_back(processorSlice);
				} catch (std::runtime_error &e) {
					err_printf(__func__, __LINE__, "Thread %s: %s\n", thread->threadDescription.c_str(), e.what());
					err_printf(__func__,
										 __LINE__,
										 "ProcessorSlice[%d] %s info: offset(%d, %d, %d), global(%d, %d, %d), local(%d, %d, %d) \n",
										 kI,
										 kernel->name.c_str(),
										 CHOOSE(kernel->offset != nullptr, kernel->offset[1], 0),
										 CHOOSE(kernel->offset != nullptr && kernel->workDim > 1, kernel->offset[1], 0),
										 CHOOSE(kernel->offset != nullptr && kernel->workDim > 2, kernel->offset[2], 0),
										 kernel->globalSize[0],
										 CHOOSE(kernel->workDim > 1, kernel->globalSize[1], 0),
										 CHOOSE(kernel->workDim > 2, kernel->globalSize[2], 0),
										 kernel->localSize[0],
										 CHOOSE(kernel->workDim > 1, kernel->localSize[1], 0),
										 CHOOSE(kernel->workDim > 2, kernel->localSize[2], 0));
				}
			}
		}
//		debug_printf(__func__, "Kernel slices execution started! Num running kernels: %d\n",
//								 runningKernels.size());
#ifdef ONE_CHUNK_AT_A_TIME
		for (auto &rk : runningKernels) {
			rk->waitForKernelEvent();
			chunks->atomicExeTime += rk->getKernelEventTime(); // for parallel level calculation
			delete rk; // clear temporary slices
			// TODO update kernelSliceTracker
		}
		runningKernels.clear();
#endif
		singleChunkT = getCurrentTime() - singleChunkT;
		thread->CT += singleChunkT;
//		debug_printf(__func__, "One slice done! singleChunkT: %f\n", singleChunkT);
		thread->addExecutedKernel(k);
	}

	for (auto &rk : runningKernels) {
		rk->waitForKernelEvent();
		chunks->atomicExeTime += rk->getKernelEventTime(); // for parallel level calculation
		delete rk; // clear temporary slices
		// TODO update kernelSliceTracker
	}

	thread->activeTime += getCurrentTime() - activeTime;
	thread->isLaunched = false;
//	debug_printf(__func__, "All slices done! active time: %f, cumulate chunks time: %f\n",
//							 thread->activeTime, thread->chunks->atomicExeTime);
	pthread_exit(nullptr);
}
#undef globalSize
#undef localSize
#undef offset

typedef std::unordered_map<std::string, std::vector<DP>> TDP_MAP;
typedef std::vector<DP> CDP;

static bool PrfIsBetter(CDP &a, CDP &b) {
	float ET_a = 0;
	float ET_b = 0;
	for (auto &dp : a) ET_a += 1 / dp.Prf;
	for (auto &dp : b) ET_b += 1 / dp.Prf;
	return ET_a < ET_b;
}

static bool ECIsBetter(CDP &a, CDP &b) {
	float EC_a = 0;
	float EC_b = 0;
	for (auto &dp : a) EC_a += dp.EC;
	for (auto &dp : b) EC_b += dp.EC;
	return EC_a < EC_b;
}

static bool ProductIsBetter(CDP &a, CDP &b) {
	float p_a = 0;
	float p_b = 0;
	for (auto &dp : a) p_a += dp.EC / dp.Prf;
	for (auto &dp : b) p_b += dp.EC / dp.Prf;
	return p_a < p_b;
}

static bool aHasLowerET(DP &a, DP &b) {
	return (1 / a.Prf < 1 / b.Prf);
}

static bool aHasLowerEC(DP &a, DP &b) {
	return (a.EC < b.EC);
}

static bool aHasLowerProduct(DP &a, DP &b) {
	return (a.EC / a.Prf < b.EC / a.Prf);
}

class CA2020 : public EDCL_Strategy {
 public:
	explicit CA2020(EDCL *edcl) : EDCL_Strategy(edcl) {}

	explicit CA2020(EDCL *edcl, ResourceManager RM) : EDCL_Strategy(edcl) { this->rm = RM; }

	void setPreference(Preference P) {
		this->preference_ = P;
		switch (preference_) {
			TEST_CASE(PERFORMANCE, CDPaIsBetter = &PrfIsBetter; DPaIsBetter = &aHasLowerEC)
			TEST_CASE(ENERGY, CDPaIsBetter = &ECIsBetter; DPaIsBetter = &aHasLowerEC)
			TEST_CASE(PRODUCT, CDPaIsBetter = &ProductIsBetter; DPaIsBetter = &aHasLowerProduct)
			default: {
				CDPaIsBetter = &PrfIsBetter;
				DPaIsBetter = &aHasLowerET;
			}
		}
	}

	// Offline scheduling, scheduling preparation
	void plan(std::vector<AtomicKernelSet> &AKSs, ResourceManager RM) override {
		priorityStages = AKSs;
		rm = RM;
		planingTime = getCurrentTime();
//	genI0();

		if (readTDPFromCSVFile(TDP_FILE_DIS_NAME, TDP) == 0) {
			planingTime = getCurrentTime() - planingTime;
			return;
		}
		// if no TDP file, create one
		TDP_MAP TDP_raw;
		if (readTDPFromCSVFile(TDP_FILE_NAME, TDP_raw) != 0) writeTDPToCSVFile(TDP_FILE_NAME, TDP_raw);
		for (auto &aks : priorityStages) {
			for (auto &k : aks->kernels) {
				std::string kernelID = genKernelID(k);
				if (TDP_raw.count(kernelID) > 0) {
					genDPs(k, TDP_raw[kernelID], TDP_FILE_NAME);
				} else {
					std::vector<DP> DPs;
					genDPs(k, DPs, TDP_FILE_NAME);
					TDP_raw[kernelID] = DPs;
				}
			}
		}
		TDP = distillTDP(TDP_FILE_DIS_NAME, TDP_raw);
		planingTime = getCurrentTime() - planingTime;
	}

	void execute(ResourceManager ResourceManager) override {
		executionTime = getCurrentTime();
		std::unordered_map<int, KernelExecutionThread *> runningThreads;
		std::queue<EDKernel> runtimeAdaptationKernels;
		std::unordered_map<EDKernel, KernelExecutionThread *> nonSatisfiedThreads;
		int stageCnt = 0;
		for (AtomicKernelSet &stage: priorityStages) {
			double overhead = getCurrentTime();
			std::vector<DP> bestCDP;
			genBestStageCDP(bestCDP, stage, 3);
			overhead = getCurrentTime() - overhead;
			debug_printf("Best CDP", "Stage %d, total: %d\n", stageCnt, priorityStages.size());
			stageCnt++;
			printCDP(bestCDP);
			debug_printf("Stg_CA2020 execution overhead", "%lf\n", overhead);
			stage->atomicExeTime = getCurrentTime();
			// Init phase
			std::vector<KernelExecutionThread *> initKernelThreads;
			for (auto &k : stage->kernels) {
//				debug_printf(__func__, "Stage Kernel %s[%d]\n", k->name.c_str(), k->kernelID);
				bool isInitKernel = false;
				for (auto &dp: bestCDP) {
					if (genKernelID(k) == dp.DP_ID && !dp.bindWithKernel) {
						initKernelThreads.push_back(new KernelExecutionThread(edcl_, rm, dp, k));
						dp.bindWithKernel = true;
						isInitKernel = true;
						break;
					}
				}
				if (!isInitKernel) runtimeAdaptationKernels.push(k);
			}
			// === PATCH ====
			// todo: uniqueL/BQueues should be empty by now, but sometimes not. Some threads from previous stage didn't finish?
			for (auto &IDQp : rm->uniqueLQueues) {
				auto q = IDQp.second;
				q->finish();
			}
			rm->uniqueLQueues.clear();

			for (auto &IDQp : rm->uniqueBQueues) {
				auto q = IDQp.second;
				q->finish();
			}
			rm->uniqueBQueues.clear();
			// === PATCH END====
			// Launching initial kernels
			for (auto &thread:initKernelThreads) {
				debug_printf(__func__, "Launching initial kernel %s, \n\t\t\t\tDP:", thread->threadDescription.c_str());
				printDP(thread->dp);
				if (!startKernelExecutionThread(thread)) {
					err_printf(__func__, __LINE__, "Initial Kernel %s starts failed! DP: ",
										 thread->threadDescription.c_str());
					printDP(thread->dp);
					printStatus(runningThreads, nonSatisfiedThreads);
					exit(EXIT_FAILURE);
				} else {
					ADD_TO_TRACK_MAP(runningThreads, thread);
				}
			}
			debug_printf(__func__, "Num adapt kernels: %d, num executing kernels: %d\n",
									 runtimeAdaptationKernels.size(), runningThreads.size());
//			sleep(10); // DEBUG, to test adapt2

			double statusTimer = getCurrentTime();
			// Online adaptation
			while (!runtimeAdaptationKernels.empty() || !runningThreads.empty() || !nonSatisfiedThreads.empty()) {
				if (fmod(getCurrentTime() - statusTimer, 5) == 0) {
					printStatus(runningThreads, nonSatisfiedThreads);
				}
				// 1. if any executing application have completed
				bool adapt1 = false;
				for (auto ktp : runningThreads) {
					auto thread = ktp.second;
					if (thread->kernelCompleted()) {
						adapt1 = true;
						joinKernelExecutionThread(thread);
						debug_printf(__func__, "Kernel %s completed!\n", thread->threadDescription.c_str());
						thread->printActiveTime();
						runningThreads.erase(thread->getTrackNum());
						delete thread;
					}
				}

				// 2. if any application arrived
				while (!runtimeAdaptationKernels.empty()) {
					auto k = runtimeAdaptationKernels.front();
					debug_printf(__func__, "Adaption started, dealing with %s\n", k->name.c_str());
					// check if k has started non-satisfied and all chunks are enqueued.
					KernelExecutionThread *adaptKernelThread = nullptr;
					if (nonSatisfiedThreads.count(k) > 0) {
						adaptKernelThread = nonSatisfiedThreads[k];
						nonSatisfiedThreads.erase(k);
						haltRunningThread(adaptKernelThread);
						if (adaptKernelThread->allSlicesEnqueued()) {
							runtimeAdaptationKernels.pop();
							ADD_TO_TRACK_MAP(runningThreads, adaptKernelThread);
							debug_printf(__func__,
													 "Kernel %s has been adopted with non-satisfied resources and all chunks has been enqueued.\n",
													 k->name.c_str());
							continue;
						} else {
							joinKernelExecutionThread(adaptKernelThread);
							debug_printf(__func__, "Kernel %s has been adopted with non-satisfied resources\n", k->name.c_str());
						}
					}
					auto DPs_K = TDP[genKernelID((k))];
					DP *bestDP_ptr = nullptr;
					for (auto &dp: DPs_K) {
						if (dp.nl >= 0 && !rm->uniqueLQueues.empty()
								&& rm->getCPUCurrentFrequency(0) != rm->getHeteroCPUFrequencies()[0][dp.fl]) {
							continue;
						}
						if (dp.nb >= 0 && !rm->uniqueBQueues.empty()
								&& rm->getCPUCurrentFrequency(6) != rm->getHeteroCPUFrequencies()[1][dp.fb]) {
							continue;
						}
						if (dp.ng >= 0 && rm->getGPUCurrentFrequency() != rm->getGPUFrequencies()[dp.fg]) {
							continue;
						}
						DP stretchDP = dp;
						reCalNewDP_Runtime(dp, runningThreads, nonSatisfiedThreads);
						if (bestDP_ptr == nullptr || DPaIsBetter(stretchDP, *bestDP_ptr)) bestDP_ptr = &dp;
					}
					DP bestDP = *bestDP_ptr;
					reCalNewDP_Runtime(bestDP, runningThreads, nonSatisfiedThreads);
					debug_printf(__func__, "Best DP for kernel %s: ", k->name.c_str());
					printDP(bestDP);
					if (adaptKernelThread == nullptr) adaptKernelThread = new KernelExecutionThread(edcl_, rm, bestDP, k);
					else adaptKernelThread->setDP(bestDP);
					if (startKernelExecutionThread(adaptKernelThread)) { // try start
						runtimeAdaptationKernels.pop();
						ADD_TO_TRACK_MAP(runningThreads, adaptKernelThread);
						debug_printf(__func__, "Thread %s starts execution with bestDP\n",
												 adaptKernelThread->threadDescription.c_str());
						break;
					} else {  // adapt2
						debug_printf(__func__, "adapt2 starts...\n");
						int freeCoresL;
						int freeCoresB;
						getNumFreeCores(freeCoresL, freeCoresB);
						debug_printf(__func__, "[%s]Free L: %d, B: %d\n", bestDP.DP_ID.c_str(), freeCoresL, freeCoresB);
						int requiredReleasedCoreL = CHOOSE(bestDP.nl >= 0, pow2(bestDP.nl) - freeCoresL, 0);
						int requiredReleasedCoreB = CHOOSE(bestDP.nb >= 0, pow2(bestDP.nb) - freeCoresB, 0);
						debug_printf(__func__, "[%s]required L: %d, B: %d\n",
												 bestDP.DP_ID.c_str(), requiredReleasedCoreL, requiredReleasedCoreB);
						std::vector<std::pair<DP *, KernelExecutionThread *>> coreToBeReleasedThreads;
						for (auto &ea : runningThreads) {
							auto runningThread = ea.second;
							if (runningThread->allSlicesEnqueued()) continue;
							auto executing_dp = &runningThread->dp;
							auto executing_kernel = runningThread->kernel;
							int left_nl = -2; // num of cores after release
							int left_nb = -2;
							if (requiredReleasedCoreL > 0 && executing_dp->nl >= 0
									&& requiredReleasedCoreL <= pow2(executing_dp->nl)) {
								left_nl = CHOOSE(requiredReleasedCoreL == pow2(executing_dp->nl),
																 -1,
																 floor(log2(pow2(executing_dp->nl) - requiredReleasedCoreL)));
							}
							if (requiredReleasedCoreB > 0 && executing_dp->nb >= 0
									&& requiredReleasedCoreB <= pow2(executing_dp->nb)) {
								left_nb = CHOOSE(requiredReleasedCoreB == pow2(executing_dp->nb),
																 -1,
																 floor(log2(pow2(executing_dp->nb) - requiredReleasedCoreB)));
							}
							if (left_nb < -1 && left_nl < -1) continue; // no enough either L or B
							if (left_nl == -1 && executing_dp->nb < 0 && executing_dp->ng < 0) continue; // only use L
							if (left_nb == -1 && executing_dp->nl < 0 && executing_dp->ng < 0) continue; // only use B
							debug_printf(__func__, "ea[%s] dp: ", executing_kernel->name.c_str());
							printDP(*executing_dp);
							debug_printf(__func__, "ea[%s](after core release) left_nl: %d, left_nb: %d\n",
													 executing_kernel->name.c_str(), left_nl, left_nb);
							if (left_nl > 2) {
								err_printf(__func__,
													 __LINE__,
													 "left_nl = %d! requiredReleasedCoreL: %d, executing_dp->nl: %d, pow2(executing_dp->nl): %d\n",
													 left_nl, requiredReleasedCoreL,
													 executing_dp->nl, pow2(executing_dp->nl));
							}
							if (left_nb > 2) {
								err_printf(__func__,
													 __LINE__,
													 "left_nb = %d! requiredReleasedCoreB: %d, executing_dp->nl: %d, pow2(executing_dp->nb): %d\n",
													 left_nb, requiredReleasedCoreB,
													 executing_dp->nb, pow2(executing_dp->nb));
							}
							assert(left_nl <= 2);
							assert(left_nb <= 2);
							float requirement;
							float tolerance = 1.2;
							switch (preference_) {
								TEST_CASE(PERFORMANCE, requirement = (1 / executing_dp->Prf) * tolerance;
										debug_printf(__func__, "ET requirement: %f\n", requirement);)
								TEST_CASE(ENERGY, requirement = executing_dp->EC * tolerance;
										debug_printf(__func__, "Energy requirement: %f\n", requirement);)
								TEST_CASE(PRODUCT, requirement = (executing_dp->EC / executing_dp->Prf) * tolerance;
										debug_printf(__func__, "Product requirement: %f\n", requirement);)
							}
							std::vector<DP *> candidateNewDPForEA;
							for (auto &dp : TDP[genKernelID(executing_kernel)]) {
								if (dp.nl >= 0 && !rm->uniqueLQueues.empty()
										&& rm->getCPUCurrentFrequency(0) != rm->getHeteroCPUFrequencies()[0][dp.fl]) { continue; }
								if (dp.nb >= 0 && !rm->uniqueBQueues.empty()
										&& rm->getCPUCurrentFrequency(6) != rm->getHeteroCPUFrequencies()[1][dp.fb]) { continue; }
								if (dp.ng >= 0 && rm->getGPUCurrentFrequency() != rm->getGPUFrequencies()[dp.fg]) { continue; }
								if (left_nl == dp.nl && dp.nb == executing_dp->nb && dp.ng == executing_dp->ng) // release L
									candidateNewDPForEA.push_back(&dp);
								if (left_nb == dp.nb && dp.nl == executing_dp->nl && dp.ng == executing_dp->ng) // release B
									candidateNewDPForEA.push_back(&dp);
								if (left_nb == dp.nb && left_nl == dp.nl && dp.ng == executing_dp->ng)
									candidateNewDPForEA.push_back(&dp);
							}
							if (candidateNewDPForEA.empty()) continue; // if no candidate dp
							DP *new_dp = nullptr;
							float dp_capability;
							for (auto dp: candidateNewDPForEA) {
								switch (preference_) {
									TEST_CASE(PERFORMANCE, dp_capability = 1 / dp->Prf)
									TEST_CASE(ENERGY, dp_capability = dp->EC)
									TEST_CASE(PRODUCT, dp_capability = dp->EC / dp->Prf)
								}
								if (dp_capability > requirement) continue;
								DP stretchedDP = *dp;
								reCalReplaceDP_Runtime(runningThread, stretchedDP, runningThreads, nonSatisfiedThreads);
								if (new_dp == nullptr || DPaIsBetter(stretchedDP, *new_dp)) new_dp = dp;
							}
							if (new_dp == nullptr) continue;
							coreToBeReleasedThreads.emplace_back(new_dp, runningThread);
							int releasedL = MAX(CHOOSE(executing_dp->nl >= 0, pow2(executing_dp->nl), 0)
																			- CHOOSE(new_dp->nl >= 0, pow2(new_dp->nl), 0), 0);
							int releasedB = MAX(CHOOSE(executing_dp->nb >= 0, pow2(executing_dp->nb), 0)
																			- CHOOSE(new_dp->nb >= 0, pow2(new_dp->nb), 0), 0);
							requiredReleasedCoreL -= releasedL;
							requiredReleasedCoreB -= releasedB;
							if (requiredReleasedCoreB == 0 && requiredReleasedCoreL == 0) break;
						}

						// stop to_be_released_threads and restart with core_released_dp
						for (auto dptp : coreToBeReleasedThreads) {
							auto runningThread = dptp.second;
							auto new_dp = *dptp.first;
							reCalReplaceDP_Runtime(runningThread, new_dp, runningThreads, nonSatisfiedThreads);
							auto old_dp = runningThread->dp;
							debug_printf(__func__, "old_dp: ");
							printDP(old_dp);
							debug_printf(__func__, "new_dp: ");
							printDP(new_dp);
							haltRunningThread(runningThread);
							joinKernelExecutionThread(runningThread);
							runningThread->setDP(new_dp);
							if (!startKernelExecutionThread(runningThread)) {
								err_printf(__func__, __LINE__, "Restarting thread failed! %s old_dp:",
													 runningThread->threadDescription.c_str());
								printDP(old_dp);
								err_printf(__func__, __LINE__, "New_dp:");
								printDP(new_dp);
								printStatus(runningThreads, nonSatisfiedThreads);
								exit(EXIT_FAILURE);
							}
							debug_printf(__func__, "adapt 2, thread %s restarted with core-released dp: ",
													 runningThread->threadDescription.c_str());
							printDP(new_dp);
						}
						// launch k
						requiredReleasedCoreL = CHOOSE(bestDP.nl >= 0, pow2(bestDP.nl), 0);
						requiredReleasedCoreB = CHOOSE(bestDP.nb >= 0, pow2(bestDP.nb), 0);
						getNumFreeCores(freeCoresL, freeCoresB);
						if (freeCoresL >= requiredReleasedCoreL && freeCoresB >= requiredReleasedCoreB) {
							adaptKernelThread->setDP(bestDP);
							runtimeAdaptationKernels.pop();
							if (!startKernelExecutionThread(adaptKernelThread)) {
								err_printf(__func__, __LINE__, "Kernel %s adapt2 failed! DP: ",
													 adaptKernelThread->threadDescription.c_str());
								printDP(adaptKernelThread->dp);
								printStatus(runningThreads, nonSatisfiedThreads);
								exit(EXIT_FAILURE);
							}
							ADD_TO_TRACK_MAP(runningThreads, adaptKernelThread);
						} else {
							int nl = CHOOSE(freeCoresL > 0, floor(log2(4 - rm->getUniqQueueUsedLittleCores())), -1);
							int nb = CHOOSE(freeCoresB > 0, floor(log2(4 - rm->getUniqQueueUsedBigCores())), -1);
							debug_printf(__func__, "CPU usage L: %d, B: %d, looking for nl: %d, nb: %d\n",
													 rm->getUniqQueueUsedLittleCores(), rm->getUniqQueueUsedBigCores(), nl, nb);
							DP *candidate_DP_ptr = nullptr;
							for (auto &dp : DPs_K) {
								if (dp.nl >= 0 && !rm->uniqueLQueues.empty()
										&& rm->getCPUCurrentFrequency(0) != rm->getHeteroCPUFrequencies()[0][dp.fl]) {
									continue;
								}
								if (dp.nb >= 0 && !rm->uniqueBQueues.empty()
										&& rm->getCPUCurrentFrequency(6) != rm->getHeteroCPUFrequencies()[1][dp.fb]) {
									continue;
								}
								if (dp.ng >= 0 && rm->getGPUCurrentFrequency() != rm->getGPUFrequencies()[dp.fg]) { continue; }
								if ((nl == dp.nl && nb == dp.nb) || (nl == dp.nl && dp.nb == -1) || (nb == dp.nb && dp.nl == -1)
										|| (dp.nl == -1 && dp.nb == -1)) {
									DP stretchedDP = dp;
									reCalNewDP_Runtime(stretchedDP, runningThreads, nonSatisfiedThreads);
									if (candidate_DP_ptr == nullptr || DPaIsBetter(stretchedDP, *candidate_DP_ptr))
										candidate_DP_ptr = &dp;
								}
							}
							DP candidate_DP = *candidate_DP_ptr;
							reCalNewDP_Runtime(candidate_DP, runningThreads, nonSatisfiedThreads);
							adaptKernelThread->setDP(candidate_DP);
							if (!startKernelExecutionThread(adaptKernelThread)) {
								err_printf(__func__, __LINE__, "Kernel %s adapt2 failed! DP: ",
													 adaptKernelThread->threadDescription.c_str());
								printDP(adaptKernelThread->dp);
								printStatus(runningThreads, nonSatisfiedThreads);
								exit(EXIT_FAILURE);
							}
							nonSatisfiedThreads[k] = adaptKernelThread;
						}
						debug_printf(__func__, "adapt2 Kernel %s started with dp:",
												 adaptKernelThread->threadDescription.c_str());
						printDP(adaptKernelThread->dp);
						debug_printf(__func__, "adapt2 end.\n");
						break;
					}
				}

				if (adapt1) {
					if (runningThreads.empty() && runtimeAdaptationKernels.empty()) continue;
					int numFreeL, numFreeB;
					getNumFreeCores(numFreeL, numFreeB);
					int free_nl = CHOOSE(numFreeL > 0, floor(log2(numFreeL)), -1);
					int free_nb = CHOOSE(numFreeB > 0, floor(log2(numFreeB)), -1);
					debug_printf(__func__, "adapt 1 starts... numFreeL: %d, numFreeB: %d. free_nl: %d, free_nb: %d \n",
											 numFreeL, numFreeB, free_nl, free_nb);
					if (numFreeL > 0 || numFreeB > 0) {
						for (auto ktp : runningThreads) {
							// Find RC, RT, REC
							auto thread = ktp.second;
							debug_printf(__func__, "Checking if thread %s can use free cores...\n",
													 thread->threadDescription.c_str());
							auto currentDP = &thread->dp;
							auto totalNC = thread->chunks->kernels.size();
							auto executedNC = thread->executedChunks.size();
							auto RC = totalNC - executedNC;
							float ATOC = (float)thread->CT / executedNC; // one chunk time
							float RT = RC * ATOC;
							float REC = (thread->dp.EC / totalNC) * RC;
							DP *bestDP_ptr = nullptr;
							DP bestDP;
							for (auto &dp : TDP[genKernelID(thread->kernel)]) {
								if (dp.nl >= 0 && !rm->uniqueLQueues.empty()
										&& rm->getCPUCurrentFrequency(0) != rm->getHeteroCPUFrequencies()[0][dp.fl]) {
									continue;
								}
								if (dp.nb >= 0 && !rm->uniqueBQueues.empty()
										&& rm->getCPUCurrentFrequency(6) != rm->getHeteroCPUFrequencies()[1][dp.fb]) {
									continue;
								}
								if (dp.ng >= 0 && rm->getGPUCurrentFrequency() != rm->getGPUFrequencies()[dp.fg]) continue;

								if (dp.nl <= free_nl && dp.nb <= free_nb) {
									DP stretchedDP = dp;
									reCalReplaceDP_Runtime(thread, stretchedDP, runningThreads, nonSatisfiedThreads);
									float RT_ = RC * ((1 / stretchedDP.Prf) / (float)stretchedDP.NC);
									float REC_ = RC * (stretchedDP.EC / (float)stretchedDP.NC);
									bool isImproved = false;
									switch (preference_) {
										TEST_CASE(PERFORMANCE, if (RT_ < RT) isImproved = true;)
										TEST_CASE(ENERGY, if (REC < REC_) isImproved = true;)
										TEST_CASE(PRODUCT, if (RT_ * REC_ < RT * REC) isImproved = true;)
									}
									if (!thread->allSlicesEnqueued() && isImproved) {
										haltRunningThread(thread);
										joinKernelExecutionThread(thread);
										DP oldDP = thread->dp;
										thread->setDP(stretchedDP);
										if (!startKernelExecutionThread(thread)) {
											thread->setDP(oldDP);
											startKernelExecutionThread(thread);
											debug_printf(__func__, "adapt 1, thread %s starts with newDP failed, restart with oldDP\n",
																	 thread->threadDescription.c_str());
											debug_printf(__func__, "numFreeL: %d, numFreeB: %d. free_nl: %d, free_nb: %d \n",
																	 numFreeL, numFreeB, free_nl, free_nb);
											continue; // try next dp
										} else {
											numFreeL -= CHOOSE(stretchedDP.nl >= 0, (int)pow2(stretchedDP.nl), 0);
											numFreeB -= CHOOSE(stretchedDP.nb >= 0, (int)pow2(stretchedDP.nb), 0);
											free_nl = CHOOSE(numFreeL > 0, floor(log2(numFreeL)), -1);
											free_nb = CHOOSE(numFreeB > 0, floor(log2(numFreeB)), -1);
											debug_printf(__func__, "adapt 1, thread %s starts with newDP success!\n",
																	 thread->threadDescription.c_str());
											printDP(thread->dp);
											debug_printf(__func__, "numFreeL: %d, numFreeB: %d. free_nl: %d, free_nb: %d \n",
																	 numFreeL, numFreeB, free_nl, free_nb);
											break;
										}
									}
								}
							}
							if (numFreeL == 0 && numFreeB == 0) break;
						}
					}
				}
			}
			// Stage finished
			stage->atomicExeTime = getCurrentTime() - stage->atomicExeTime;
		}
		executionTime = getCurrentTime() - executionTime;
	}

// Releasing resources, free up tmp memory
	void finish(ResourceManager ResourceManager) override {

	}

 protected:
	std::vector<DP> I0_DPV;
	std::string I0_id = "I0";
	TDP_MAP TDP;
	ResourceManager rm{};
	Preference preference_ = PERFORMANCE;
	bool (*CDPaIsBetter)(CDP &a, CDP &b) = &PrfIsBetter;
	bool (*DPaIsBetter)(DP &a, DP &b) = &aHasLowerEC;

	void printStatus(std::unordered_map<int, KernelExecutionThread *> &runningThreads,
									 std::unordered_map<EDKernel, KernelExecutionThread *> &nonSatisfiedThreads) {
		debug_printf("[Running Status]", "\n");
		debug_printf("===Running Status===", "CPU usage L: %d, B: %d\n",
								 rm->getUniqQueueUsedLittleCores(), rm->getUniqQueueUsedBigCores());
		debug_printf("===Running Status===", "Running threads: \n");
		for (auto &ktp:runningThreads) {
			auto rth = ktp.second;
			debug_printf("===Running Status===",
									 " Thread %s, total chunks: %d, executed chunks: %d, chunk time: %f, \n\t\t\t\tDP: ",
									 rth->threadDescription.c_str(), rth->chunks->kernels.size(),
									 rth->getNumExecutedKernels(), rth->CT);
			printDP(rth->dp);
		}
		for (auto &ktp:nonSatisfiedThreads) {
			auto rth = ktp.second;
			debug_printf("===Running Status===",
									 " Thread %s, total chunks: %d, executed chunks: %d, chunk time: %f, \n\t\t\t\tDP: ",
									 rth->threadDescription.c_str(), rth->chunks->kernels.size(),
									 rth->getNumExecutedKernels(), rth->CT);
			printDP(rth->dp);
		}
	}

	void getNumFreeCores(int &freeCoresL,
											 int &freeCoresB) {
		freeCoresL = 4 - rm->getUniqQueueUsedLittleCores();
		freeCoresB = 4 - rm->getUniqQueueUsedBigCores();
	}

	static bool startKernelExecutionThread(KernelExecutionThread *thread) {
//		debug_printf(__func__, "launching %s...\n", thread->threadDescription.c_str());
		if (thread->isLaunched) return true;
		if (!thread->prepareExecution()) return false;
		pthread_mutex_trylock(&thread->runningMutex); // NOLINT(bugprone-unused-return-value)
		int rc;
		if ((rc = pthread_create(&thread->pthread, nullptr, KernelExecutionThreadFunction, thread))) {
			fprintf(stderr, "[%s]error: pthread_create, rc: %d\n", __func__, rc);
			return false;
		} else thread->isLaunched = true;
//		debug_printf(__func__, "%s launched!\n", thread->threadDescription.c_str());
		return true;
	}

	static void haltRunningThread(KernelExecutionThread *thread) {
		pthread_mutex_unlock(&thread->runningMutex);
	}

	static void joinKernelExecutionThread(KernelExecutionThread *executionThread) {
		if (!executionThread->isLaunched) return;
		pthread_join(executionThread->pthread, nullptr);
		executionThread->releaseQueues();
		executionThread->isLaunched = false;
	}

	static float getStretchedGPUTime(std::unordered_map<int, KernelExecutionThread *> &runningThreads,
																	 std::unordered_map<EDKernel, KernelExecutionThread *> &nonSatisfiedThreads,
																	 KernelExecutionThread *excludeThread) {
		float runningGPUT = 0;
		for (auto ktp : runningThreads) {
			auto thread = ktp.second;
			if (thread->kernelCompleted()) continue;
			if (excludeThread != nullptr && thread->getTrackNum() == excludeThread->getTrackNum()) continue;
			if (thread->dp.ng >= 0) runningGPUT += 1 / thread->dp.Prf;
		}
		for (auto ktp : nonSatisfiedThreads) {
			auto thread = ktp.second;
			if (thread->kernelCompleted()) continue;
			if (excludeThread != nullptr && thread->getTrackNum() == excludeThread->getTrackNum()) continue;
			if (thread->dp.ng >= 0) runningGPUT += 1 / thread->dp.Prf;
		}
		return runningGPUT;
	}

	void reCalNewDP_Runtime(DP &dp,
													std::unordered_map<int, KernelExecutionThread *> &runningThreads,
													std::unordered_map<EDKernel, KernelExecutionThread *> &nonSatisfiedThreads) {
		float gSET = getStretchedGPUTime(runningThreads, nonSatisfiedThreads, nullptr);
		reCalDP(dp, gSET);
	}

	void reCalReplaceDP_Runtime(KernelExecutionThread *thread,
															DP &dp,
															std::unordered_map<int, KernelExecutionThread *> &runningThreads,
															std::unordered_map<EDKernel, KernelExecutionThread *> &nonSatisfiedThreads) {
		float gSET = getStretchedGPUTime(runningThreads, nonSatisfiedThreads, thread);
		reCalDP(dp, gSET);
	}

	void reCalDP(DP &dp, float gSET) {
		DP dp_org = dp;
		assert(gSET >= 0);
		if (dp.ng >= 0) gSET += 1 / dp.Prf;
		std::vector<float> prfs{0, 0, 0};
		prfs[0] = CHOOSE(dp.nl >= 0, dp.Prf, 0);
		prfs[1] = CHOOSE(dp.nb >= 0, dp.Prf, 0);
		prfs[2] = CHOOSE(dp.ng >= 0, 1 / gSET, 0);
		float sumPrfs = 0;
		size_t totalCK = 0;
		for (int kI = 0; kI < 3; ++kI) {
			sumPrfs += prfs[kI];
			totalCK += dp.CK[kI];
		}
		size_t recalCK = 0;
		int maxCK_device;
		size_t maxCK = 0;
		std::vector<size_t> oldCK = dp.CK;
		for (int i = 0; i < 3; ++i) {
			dp.CK_ratio[i] = prfs[i] / sumPrfs;
			dp.CK[i] = ceil(totalCK * dp.CK_ratio[i]);
			recalCK += dp.CK[i];
			if (dp.CK[i] > maxCK) {
				maxCK = dp.CK[i];
				maxCK_device = i;
			}
		}
//	  flushed_printf("\n");
//	  debug_printf(__func__, "totalCK: %d, recalCK: %d\n", totalCK, recalCK);
		float oldET = 1 / dp.Prf;
		if (recalCK != totalCK) { dp.CK[maxCK_device] += totalCK - recalCK; }
//	  debug_printf(__func__, "Final CK(old, new): ");
//	  for (int i = 0; i < 3; ++i) {
//		flushed_printf("%d(%d, %d) ", i, oldCK[i], dp.CK[i]);
//	  }
//	  flushed_printf("\n");

		// recal ET
		float newOldCKRatioL = CHOOSE(oldCK[0] > 0, (float)dp.CK[0] / (float)oldCK[0], 0.0);
		float newOldCKRatioB = CHOOSE(oldCK[1] > 0, (float)dp.CK[1] / (float)oldCK[1], 0.0);
		float newOldCKRatioG = CHOOSE(oldCK[2] > 0, (float)dp.CK[2] / (float)oldCK[2], 0.0);
		float newET_L = oldET * newOldCKRatioL;
		float newET_B = oldET * newOldCKRatioB;
		float newET_G = oldET * newOldCKRatioG;
		float newET = MAX(newET_L, MAX(newET_B, newET_G));
		dp.Prf = 1 / newET;
		if (newET <= 0) {
			debug_printf(__func__, "prfs[0]=%f, prfs[1]=%f, prfs[2]=%f\n", prfs[0], prfs[1], prfs[2]);
			debug_printf(__func__, "gSET: %f, newOldCKRatioL: %f, newOldCKRatioB: %f, newOldCKRatioG: %f\n",
									 gSET, newOldCKRatioL, newOldCKRatioB, newOldCKRatioG);
			debug_printf(__func__, "newET_L: %f, newET_B: %f, newET_G: %f, newET: %f, dp.Prf: %f\n",
									 newET_L, newET_B, newET_G, newET, dp.Prf);
			printDP(dp_org);
			printDP(dp);
			assert(newET > 0);
		}

		// recal EC
//	  float i0 = CHOOSE(dp.nl >= 0, I0[dp.fl], I0[3]);
//	  float lambdaL0 = CHOOSE(dp.nl >= 0, LAMBDA_L[dp.fl][0], LAMBDA_L[3][0]);
//	  float f0 = CHOOSE(dp.nl >= 0,
//						(float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f,
//						(float)rm->getHeteroCPUFrequencies()[0][3] / 1000000.0f);
		float fl = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f, 0);
		float fb = CHOOSE(dp.nb >= 0, (float)(rm->getHeteroCPUFrequencies()[1].at(dp.fb)) / 1000000.0f, 0);
		float fg = CHOOSE(dp.ng >= 0, (float)(rm->getGPUFrequencies().at(dp.fg)) / 1000.0f, 0);
		float lambdaL = 0;
		if (dp.nl >= 0) lambdaL = LAMBDA_L[dp.fl][dp.nl];
		float lambdaB = 0;
		if (dp.nb >= 0) lambdaB = LAMBDA_B[dp.fb][dp.nb];
		float lambdaG = 0;
		if (dp.ng >= 0) lambdaG = LAMBDA_G[dp.fg];
		float I_sqr = (dp.EC / oldET) / (lambdaL * fl + lambdaB * fb + lambdaG * fg);
		dp.EC = I_sqr *
				(newOldCKRatioL * lambdaL * fl + newOldCKRatioB * lambdaB * fb + newOldCKRatioG * lambdaG * fg) * newET;
//	  if (dp.EC < 0) {
//		debug_printf(__func__,
//					 "EC: %f, I_sqr: %f, newOldCKRatioL: %f, lambdaL: %f, fl: %f, newET_L: %f, newOldCKRatioB: %f, lambdaB: %f, fb: %f, newET_B: %f, newOldCKRatioG: %f, lambdaG: %f, fg: %f, newET_G: %f, i0: %f, lambdaL0: %f, f0: %f, newET: %f\n",
//					 dp.EC,
//					 I_sqr,
//					 newOldCKRatioL,
//					 lambdaL,
//					 fl,
//					 newET_L,
//					 newOldCKRatioB,
//					 lambdaB,
//					 fb,
//					 newET_B,
//					 newOldCKRatioG,
//					 lambdaG,
//					 fg,
//					 newET_G,
//					 i0,
//					 lambdaL0,
//					 newET);
//	  }

	}

	void reCalCDP(std::vector<DP> &CDP) {
		float gSET = 0;
		for (auto &dp_cdp: CDP) {
			if (dp_cdp.ng >= 0) gSET += 1 / dp_cdp.Prf;
//	  debug_printf(__func__, "dp_cdp ET: %f\n", dp_cdp.Prf);
		}
//	debug_printf(__func__, "gSET: %f\n", gSET);

		for (auto &dp_cdp: CDP) {
			reCalDP(dp_cdp, gSET);
		}
	}

	void genBestStageCDP(std::vector<DP> &CDP,
											 AtomicKernelSet &stage, int minNumDP) {
		CDP.clear();
		std::vector<std::vector<DP> *> stageDPVs;
		for (int i = 0; i < stage->kernels.size() && i < minNumDP; i++) {
			stageDPVs.push_back(&TDP[genKernelID(stage->kernels[i])]);
		}
		int numK = stageDPVs.size();
		std::vector<int> counter(numK);
		std::vector<int> counterLim(numK);
//	debug_printf(__func__, "counterLim[");
		for (int i = 0; i < numK; i++) {
			counter[i] = 0;
			counterLim[i] = stageDPVs[i]->size();
//	  flushed_printf(" %d ", counterLim[i]);
		}
//	flushed_printf("]\n");

		bool allDPsAccessed = false;
		while (!allDPsAccessed) {
			++counter[0];
			for (int i = 1; i < numK; ++i) {
				if (counter[i - 1] >= counterLim[i - 1]) {
					counter[i - 1] = 0;
					counter[i]++;
				}
			}
			allDPsAccessed = true;
			for (int i = 0; i < numK; ++i) {
				if (counter[i] < counterLim[i] - 1) {
					allDPsAccessed = false;
					break;
				}
			}

			int fL = -1;
			int fB = -1;
			int fG = -1;
			bool skip = false;
			for (int i = 0; i < numK; ++i) {
				DP *dp_p = &stageDPVs[i]->at(counter[i]);
				fL = CHOOSE(dp_p->nl >= 0 && fL == -1, fL = dp_p->fl, fL);
				fB = CHOOSE(dp_p->nb >= 0 && fB == -1, fB = dp_p->fb, fB);
				fG = CHOOSE(dp_p->ng >= 0 && fG == -1, fG = dp_p->fg, fG);
				skip |= dp_p->nl >= 0 && dp_p->fl != fL;
				skip |= dp_p->nb >= 0 && dp_p->fb != fB;
				skip |= dp_p->ng >= 0 && dp_p->fg != fG;
				if (skip) break;
			}
			if (skip) continue;

			std::vector<DP> CDP_candidate;
			int numL = 0;
			int numB = 0;
			for (int i = 0; i < numK; i++) {
				DP *dp_p = &stageDPVs[i]->at(counter[i]);
				int dp_numL = CHOOSE(dp_p->nl >= 0, pow2(dp_p->nl), 0);
				int dp_numB = CHOOSE(dp_p->nb >= 0, pow2(dp_p->nb), 0);
				if (dp_numL + numL > 4) continue;
				if (dp_numB + numB > 4) continue;
				numL = dp_numL + numL;
				numB = dp_numB + numB;
				CDP_candidate.push_back(*dp_p);
			}
			if (CDP_candidate.size() >= minNumDP || stage->kernels.size() == CDP_candidate.size()) {
				reCalCDP(CDP_candidate);
				if (CDP.empty() || CDPaIsBetter(CDP_candidate, CDP)) { CDP = CDP_candidate; }
			}
		}
	}

	void genAllStageCDP(std::vector<std::vector<DP> > &AllCDPs,
											AtomicKernelSet &stage, int minNumDP) {
		std::vector<std::vector<DP> *> stageDPVs;
		for (auto &k : stage->kernels) {
			stageDPVs.push_back(&TDP[genKernelID(k)]);
		}
		int numK = stageDPVs.size();
		std::vector<int> counter(numK);
		std::vector<int> counterLim(numK);
		debug_printf(__func__, "counterLim[");
		for (int i = 0; i < numK; i++) {
			counter[i] = 0;
			counterLim[i] = stageDPVs[i]->size();
			flushed_printf(" %d ", counterLim[i]);
		}
		flushed_printf("]\n");

		bool allDPsAccessed = false;
		while (!allDPsAccessed) {
			++counter[0];
			for (int i = 1; i < numK; ++i) {
				if (counter[i - 1] >= counterLim[i - 1]) {
					counter[i - 1] = 0;
					counter[i]++;
				}
			}
			allDPsAccessed = true;
			for (int i = 0; i < numK; ++i) {
				if (counter[i] < counterLim[i] - 1) {
					allDPsAccessed = false;
					break;
				}
			}

			int fL = -1;
			int fB = -1;
			int fG = -1;
			bool skip = false;
			for (int i = 0; i < numK; ++i) {
				DP *dp_p = &stageDPVs[i]->at(counter[i]);
				fL = CHOOSE(dp_p->nl >= 0 && fL == -1, fL = dp_p->fl, fL);
				fB = CHOOSE(dp_p->nb >= 0 && fB == -1, fB = dp_p->fb, fB);
				fG = CHOOSE(dp_p->ng >= 0 && fG == -1, fG = dp_p->fg, fG);
				skip |= dp_p->nl >= 0 && dp_p->fl != fL;
				skip |= dp_p->nb >= 0 && dp_p->fb != fB;
				skip |= dp_p->ng >= 0 && dp_p->fg != fG;
				if (skip) break;
			}
			if (skip) continue;

			std::vector<DP> CDP_candidate;
			int numL = 0;
			int numB = 0;
			for (int i = 0; i < numK; i++) {
				DP *dp_p = &stageDPVs[i]->at(counter[i]);
				int dp_numL = CHOOSE(dp_p->nl >= 0, pow2(dp_p->nl), 0);
				int dp_numB = CHOOSE(dp_p->nb >= 0, pow2(dp_p->nb), 0);
				if (dp_numL + numL > 4) continue;
				if (dp_numB + numB > 4) continue;
				numL = dp_numL + numL;
				numB = dp_numB + numB;
				CDP_candidate.push_back(*dp_p);
			}

			if (CDP_candidate.size() >= minNumDP || stage->kernels.size() == CDP_candidate.size()) {
				reCalCDP(CDP_candidate);
				AllCDPs.push_back(CDP_candidate);
			}
		}
	}

	void genI0() {
		TDP_MAP I0_TDP;
		if (readTDPFromCSVFile("I0.csv", I0_TDP) == 0) {
			I0_DPV = I0_TDP[I0_id];
			return;
		}
		auto freq_L = rm->getHeteroCPUFrequencies()[0];
		auto freq_B = rm->getHeteroCPUFrequencies()[1];
		auto freq_G = rm->getGPUFrequencies();
		int t_nfL = freq_L.size(); // total num little core frequency
		int t_nfB = freq_B.size();
		int t_nfG = freq_G.size();
		int sleepTime = 120;
		float avgCurrent = 0;
		std::vector<float> current;
		// sample a value of current every 100 usec as defined in EDCL_ResourceManager.h
		current.reserve((sleepTime * 1000000) / 100);
		for (int nfL = 0; nfL < t_nfL; ++nfL) {
			for (int nfB = 0; nfB < t_nfB; ++nfB) {
				for (int nfG = 0; nfG < t_nfG; ++nfG) {
					DP dp;
					dp.DP_ID = I0_id;
					dp.Prf = pow(sleepTime, -1);
					dp.fb = nfB;
					dp.fl = nfL;
					dp.fg = nfG;
					rm->setCPUFrequency(0, dp.fl);
					rm->setCPUFrequency(rm->getNumCPU() - 1, dp.fb);
					rm->setGPUFrequency(dp.fg);
					sleep(sleepTime);
					if (rm->startCurrentRecording(&current) == 0) {
						sleep(sleepTime);
						rm->stopCurrentRecording();
					}
					debug_printf(__func__, "Num records: %d. ", current.size());
					for (auto &c:current) {
						avgCurrent += c;
					}
					avgCurrent /= current.size();
					dp.EC = avgCurrent;
					PRINT_CYAN
					flushed_printf("Freq: [%.4f, %.4f, %.4f]. ",
												 freq_L[nfL] / 1000000.0,
												 freq_B[nfB] / 1000000.0,
												 freq_G[nfG] / 1000.0);
					flushed_printf("Avg current: %f\n", avgCurrent);
					PRINT_DEFAULT
					current.clear();
					avgCurrent = 0;
					I0_DPV.push_back(dp);
				}
			}
		}
		I0_TDP[I0_id] = I0_DPV;
		writeTDPToCSVFile("I0.csv", I0_TDP);
	}

	static void getPerf(DP &dp, std::vector<float> &Prfs, std::vector<DP> &DPs) {
		// single device case
		if (dp.nl >= 0 && dp.nb < 0 && dp.ng < 0) {
			Prfs[0] = 1;
			return;
		} else if (dp.nl < 0 && dp.nb >= 0 && dp.ng < 0) {
			Prfs[1] = 1;
			return;
		} else if (dp.nl < 0 && dp.nb < 0 && dp.ng >= 0) {
			Prfs[2] = 1;
			return;
		}

		// general case
		for (auto &dataPoint: DPs) { // find the performance indicator of the kernel executing on the only device
			if (dp.nl >= 0 && dp.nl == dataPoint.nl && dp.fl == dataPoint.fl
					&& dataPoint.nb < 0 && dataPoint.ng < 0) {
				Prfs[0] = dataPoint.Prf;
			}
			if (dp.nb >= 0 && dp.nb == dataPoint.nb && dp.fb == dataPoint.fb
					&& dataPoint.nl < 0 && dataPoint.ng < 0) {
				Prfs[1] = dataPoint.Prf;
			}
			if (dp.ng >= 0 && dp.ng == dataPoint.ng && dp.fg == dataPoint.fg
					&& dataPoint.nb < 0 && dataPoint.nl < 0) {
				Prfs[2] = dataPoint.Prf;
			}
			bool allSet = true;
			for (auto &prf : Prfs) {
				if (prf == 0) allSet = false;
			}
			if (allSet) break;
		}
	}


// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()
	static std::string genKernelID(EDKernel k) {
		std::string kernel_id = k->name;
		for (int kI = 0; kI < k->workDim; ++kI) {
			kernel_id += "_" + std::to_string(k->globalSize[kI]) + "_" + std::to_string(k->localSize[kI]);
		}
		return kernel_id;
	}

	void distributeWorkload(AtomicKernelSet &kernelSlices,
													DP &dp,
													std::vector<EDQueue> &queues,
													OCLKernelExecutionThread **packs, int heteroLevel) {
//	printDP(dp);
		// Launch kernel execution threads
		for (int kI = 0; kI < heteroLevel; ++kI) {
//	  debug_printf(__func__,
//				   "%d: %d, %d, %d\n",
//				   kI,
//				   dp.CK_ratio[kI] > 0,
//				   queues.at(kI) != nullptr,
//				   packs[kI] == nullptr);
			// use CK_ration[kI] to test which devices are used
			if (dp.CK_ratio[kI] > 0 && queues.at(kI) != nullptr && packs[kI] == nullptr) {
				packs[kI] = new OCLKernelExecutionThread(edcl_, queues[kI]);
//		debug_printf(__func__, "pack %d created!\n", kI);
//		rm->startOCLKernelExecutionThread(packs[kI]);
			} //else packs[kI] = nullptr;
		}
		for (auto &k: kernelSlices->kernels) { // for each slice
//	  debug_printf(__func__, "Kernel_in ");
//	  k->printInfo();
			// slice on the largest dim
			int largestDim = 0;
			for (int kJ = 0; kJ < k->workDim; ++kJ) {
				if (k->globalSize[kJ] / k->localSize[kJ] > k->globalSize[largestDim] / k->localSize[largestDim])
					largestDim = kJ;
			}
			size_t totalNumWG = k->globalSize[largestDim] / k->localSize[largestDim];
//	  debug_printf(__func__, "largestDim is %d, totalNumWG: %d\n", largestDim, totalNumWG);

			int lastUsedDevice = 0;
			for (int kI = 0; kI < heteroLevel; ++kI) {
				if (dp.CK_ratio[kI] > 0) lastUsedDevice = kI;
			}

			// DEBUG
			size_t distributedNumWG = 0;
			// DEBUG END

			std::unique_ptr<size_t[]> k_offset = std::make_unique<size_t[]>(k->workDim);
			std::unique_ptr<size_t[]> k_global = std::make_unique<size_t[]>(k->workDim);
			memset(k_offset.get(), 0, k->workDim * sizeof(size_t));
			memset(k_global.get(), 0, k->workDim * sizeof(size_t));
			for (int kI = 0; kI < heteroLevel; ++kI) { // kI: loop through processors
				if (dp.CK_ratio[kI] > 0 && packs[kI] != nullptr) {
					EDKernel processorSlice = nullptr;
					for (int kJ = 0; kJ < k->workDim; ++kJ) { //kJ: loop through workDim
						if (kJ == largestDim) {
							size_t numSliceWG = floor((float)totalNumWG * dp.CK_ratio[kI]);
							if (numSliceWG > 0 || (kI == lastUsedDevice && k_offset[kJ] < k->globalSize[kJ])) {
								processorSlice = edcl_->createKernelCopy(k);
								k_global[kJ] = numSliceWG * k->localSize[kJ];
							}
						} else k_global[kJ] = k->globalSize[kJ];
					}

					if (processorSlice != nullptr) {
						UNIQUE_PTR_SIZE_T(sliceGlobal, processorSlice->workDim);
						UNIQUE_PTR_SIZE_T(sliceOffset, processorSlice->workDim);
						processorSlice->setGlobalSize(sliceGlobal.get());
						processorSlice->setOffset(sliceOffset.get());
						for (int kJ = 0; kJ < k->workDim; ++kJ) { // loop through work-dim
							processorSlice->offset[kJ] = k_offset[kJ] + k->getOffset()[kJ];
//			  debug_printf(__func__, "Tracking: k_offset[%d] = %d, k_global[%d] = %d\n",
//						   kJ, k_offset[kJ], kJ, k_global[kJ]);
							// if it's the last device and there're more than global size work left, the last device handles the rest
							if (kI == lastUsedDevice && k->globalSize[kJ] - k_offset[kJ] > 0)
								processorSlice->globalSize[kJ] = k->globalSize[kJ] - k_offset[kJ];
							else
								processorSlice->globalSize[kJ] = k_global[kJ];
//			  debug_printf(__func__, "%s slice offset[%d]: %d, global[%d]: %d. LargestDim: %d\n",
//						   processorSlice->name.c_str(),
//						   kJ, processorSlice->offset[kJ],
//						   kJ, processorSlice->globalSize[kJ], largestDim);
//			  debug_printf(__func__, "Device[%d] get %d WGs\n", kI, dp.CK[kI]);
							if (kJ == largestDim) {
								k_offset[kJ] += k_global[kJ];
							}
						}
						dp.CK[kI] = processorSlice->globalSize[largestDim] / processorSlice->localSize[largestDim];
						distributedNumWG += dp.CK[kI];
						edcl_->confirmExeEnv(packs[kI]->queue, processorSlice);
						packs[kI]->kernels.enqueue(processorSlice);
//			debug_printf(__func__, "Slice ");
//			processorSlice->printInfo();
					}
				}
			}
			// DEBUG
			if (distributedNumWG != totalNumWG) {
				err_printf(__func__, __LINE__, "distributedNumWG(%d) != totalNumWG(%d)! LargestDim: %d\n",
									 distributedNumWG, totalNumWG, largestDim);
				printDPVertically(dp, rm);
				EDKernel kernel = k;
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
				for (int kI = 0; kI < heteroLevel; ++kI) {
					if (packs[kI] != nullptr) {
						while (packs[kI]->kernels.try_dequeue(kernel)) {
							kernel->printInfo();
						}
					}
				}
				exit(-1);
			}
			// DEBUG END
		}
	}
#undef globalSize
#undef localSize
#undef offset

	double concurrentWorkExecution(DP &dp, OCLKernelExecutionThread **packs, int heteroLevel) {
		// Launch kernel execution threads
		double time = getCurrentTime();
		for (int kI = 0; kI < heteroLevel; ++kI) {
			if (dp.CK_ratio[kI] > 0) {
//		debug_printf(__func__, "Device[%d] %s thread started, get num kernels: %d\n",
//					 kI, packs[kI]->packDes.c_str(), packs[kI]->kernels.size_approx());
				rm->startOCLKernelExecutionThread(packs[kI]);
			}
		}
//	debug_printf(__func__, "OCLKernel exe thread started!\n");

		for (int kI = 0; kI < heteroLevel; ++kI) {
			if (packs[kI] != nullptr) {
				rm->unlockRunningThread(packs[kI]);
			}
		}

		for (int kI = 0; kI < heteroLevel; ++kI) {
			if (dp.CK_ratio[kI] > 0) {
				rm->joinOCLKernelExecutionThread(packs[kI]);
//		debug_printf(__func__, "Device[%d] %s thread joined, num executed kernels: %d\n",
//					 kI, packs[kI]->packDes.c_str(), packs[kI]->executedKernels.size());
//		debug_printf(__func__, "%s pack concurrentQ: %d\n",
//					 packs[kI]->packDes.c_str(), packs[kI]->kernels.size_approx());
			}
		}
		time = getCurrentTime() - time;
//	debug_printf(__func__, "All OCLKernel exe threads joined! Total time: %f\n", time);
		return time;
	}

	double benchmarkKernelSlices(DP &dp, AtomicKernelSet &kernelSlices, std::vector<EDQueue> &queues) {
		OCLKernelExecutionThread *packs[3] = {nullptr, nullptr, nullptr};
		distributeWorkload(kernelSlices, dp, queues, packs, 3);
		double time = concurrentWorkExecution(dp, packs, 3);
		for (auto &pack : packs) {
			if (pack != nullptr) {
				for (auto &k : pack->executedKernels) delete k;
			}
			delete pack;
		}
		return time;
	}

// once NC is determined, find out dp.prf, dp.ec
	void setDP(DP &dp, std::shared_ptr<InternalAtomicKernelSet_> kernelSlices, std::vector<EDQueue> &queues) {
//	debug_printf(__func__, "start\n");
		OCLKernelExecutionThread *packs[3] = {nullptr, nullptr, nullptr};

		distributeWorkload(kernelSlices, dp, queues, packs, 3);
		double time = concurrentWorkExecution(dp, packs, 3);

		int numLoop = MAX(BENCHMARK_KERNEL_SLICE_LOOP_NUM, ceil(SET_DP_TIME / time));
//	debug_printf(__func__, "One run time: %lf, numLoop: %d\n", time, numLoop);
		std::vector<float> current;
		current.reserve(100000);

		// store kernel slices
		std::vector<std::vector<EDKernel >> kernels;
		kernels.resize(3);
		for (int kI = 0; kI < 3; ++kI) {
			if (dp.CK_ratio[kI] > 0) {
				kernels[kI] = packs[kI]->executedKernels;
//		debug_printf(__func__, "kernels[%d] size: %d\n", kI, kernels[kI].size());
			}
		}

//	debug_printf(__func__, "Cooling down...\n");
//  DEBUG
//	sleep(SET_DP_TIME);

		for (int kI = 0; kI < 3; ++kI) {
			if (dp.CK_ratio[kI] > 0) {
				rm->startOCLKernelExecutionThread(packs[kI]);
			}
		}

//	debug_printf(__func__, "OCL executing thread started\n");
		rm->startCurrentRecording(&current);
		time = getCurrentTime();
		for (int i = 0; i < numLoop; ++i) {
			for (int kI = 0; kI < 3; ++kI) {
				if (dp.CK_ratio[kI] > 0) {
					for (EDKernel &k : kernels[kI]) {
						EDKernel k_copy = edcl_->createKernelCopy(k);
//			edcl_->confirmExeEnv(packs[kI]->queue, k_copy);
						packs[kI]->kernels.enqueue(k_copy);
					}
					packs[kI]->executedKernels.clear();
				}
			}
		}

		for (auto &pack : packs) {
			if (pack != nullptr) {
				rm->unlockRunningThread(pack);
			}
		}
//	debug_printf(__func__, "All kernels enqueued! Running lock unlocked\n");

		for (int kI = 0; kI < 3; ++kI) {
			if (dp.CK_ratio[kI] > 0) {
				rm->joinOCLKernelExecutionThread(packs[kI]);
//		debug_printf(__func__, "Device[%d] %s thread joined, num executed kernels: %d\n",
//					 kI, packs[kI]->packDes.c_str(), packs[kI]->executedKernels.size());
			}
		}

		time = getCurrentTime() - time;
		rm->stopCurrentRecording();
//	debug_printf(__func__, "Total execution time: %f, avg time: %f\n", time, time / numLoop);

		float avgC = 0;
		for (auto &c : current) {
			avgC += c;
		}
		avgC /= current.size();
		time = time / numLoop;
		dp.Prf = (float)(1.0 / time);
		calculateEC(avgC, (float)time, dp);
		for (auto &pack : packs) {
			if (pack != nullptr) {
				for (auto &k : pack->executedKernels) delete k;
			}
			delete pack;
		}

	}

	[[deprecated("Used only when converting old EC(fixed i0) to the latest")]]
	__unused void recalEC(float I, float time, DP &dp) {
//	float i0 = CHOOSE(dp.nl >= 0, I0[dp.fl], I0[3]);
//	float lambdaL0 = CHOOSE(dp.nl >= 0, LAMBDA_L[dp.fl][0], LAMBDA_L[3][0]);
//	float f0 = CHOOSE(dp.nl >= 0,
//					  (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f,
//					  (float)rm->getHeteroCPUFrequencies()[0][3] / 1000000.0f);

		float fl = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f, 0);
		float fb = CHOOSE(dp.nb >= 0, (float)(rm->getHeteroCPUFrequencies()[1].at(dp.fb)) / 1000000.0f, 0);
		float fg = CHOOSE(dp.ng >= 0, (float)(rm->getGPUFrequencies().at(dp.fg)) / 1000.0f, 0);
		float lambdaL = 0;
		if (dp.nl >= 0) lambdaL = LAMBDA_L[dp.fl][dp.nl];
		float lambdaB = 0;
		if (dp.nb >= 0) lambdaB = LAMBDA_B[dp.fb][dp.nb];
		float lambdaG = 0;
		if (dp.ng >= 0) lambdaG = LAMBDA_G[dp.fg];
//	dp.EC = (I * I * (lambdaL * fl + lambdaB * fb + lambdaG * fg)
//		- i0 * i0 * lambdaL0 * f0) * time;
		dp.EC = I * I * (lambdaL * fl + lambdaB * fb + lambdaG * fg) * time;
//	dp.EC *= 1000000;
//	if (dp.EC < 0) {
//	  debug_printf(__func__,
//				   "EC: %f, I: %f, lambdaL: %f, fl: %f, lambdaB: %f, fb: %f, lambdaG: %f, fg: %f, i0: %f, lambdaL0: %f, f0: %f, time: %f\n",
//				   dp.EC,
//				   I,
//				   lambdaL,
//				   fl,
//				   lambdaB,
//				   fb,
//				   lambdaG,
//				   fg,
//				   i0,
//				   lambdaL0,
//				   f0,
//				   time);
//	}
	}

	void calculateEC(float I, float time, DP &dp) {
//		auto currLFreq = rm->getCPUCurrentFrequency(0);
//		int currLFreqIdx;
//		auto LFreqArr = rm->getHeteroCPUFrequencies()[0];
//		for (int i = 0; i < LFreqArr.size(); i++) {
//			if (LFreqArr[i] == currLFreq) currLFreqIdx = i;
//		}
//		float i0 = I0[currLFreqIdx];
//		float lambdaL0 = LAMBDA_L[currLFreqIdx][0];
//
//		float f0 = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f,
//											rm->getHeteroCPUFrequencies()[0].back() / 1000000.0f);
		float fl = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f, 0);
		float fb = CHOOSE(dp.nb >= 0, (float)(rm->getHeteroCPUFrequencies()[1].at(dp.fb)) / 1000000.0f, 0);
		float fg = CHOOSE(dp.ng >= 0, (float)(rm->getGPUFrequencies().at(dp.fg)) / 1000.0f, 0);
		float lambdaL = 0;
		if (dp.nl >= 0) lambdaL = LAMBDA_L[dp.fl][dp.nl];
		float lambdaB = 0;
		if (dp.nb >= 0) lambdaB = LAMBDA_B[dp.fb][dp.nb];
		float lambdaG = 0;
		if (dp.ng >= 0) lambdaG = LAMBDA_G[dp.fg];
//		dp.EC = (I * I * (lambdaL * fl + lambdaB * fb + lambdaG * fg)
//				- i0 * i0 * lambdaL0 * f0) * time;
		dp.EC = I * I * (lambdaL * fl + lambdaB * fb + lambdaG * fg) * time;
		dp.EC *= 1000;
//		if (dp.EC < 0) {
//			debug_printf(__func__,
//									 "EC: %f, I: %f, lambdaL: %f, fl: %f, lambdaB: %f, fb: %f, lambdaG: %f, fg: %f, i0: %f, lambdaL0: %f, f0: %f, time: %f\n",
//									 dp.EC,
//									 I,
//									 lambdaL,
//									 fl,
//									 lambdaB,
//									 fb,
//									 lambdaG,
//									 fg,
//									 i0,
//									 lambdaL0,
//									 f0,
//									 time);
//		}
//		assert(dp.EC > 0);
	}

//	float queryI0(DP &dp) {
//		for (auto &i0: I0_DPV) {
//			if (dp.fl == i0.fl && dp.fb == i0.fb && dp.fg == i0.fg) return i0.EC;
//		}
//		return 0;
//	}

	void findNC(EDKernel kernel, DP &dp, std::vector<EDQueue> &queues) {
		uint sliceFactor = 0;
		auto aks = edcl_->slicingKernelPow2_extra_equal(kernel, sliceFactor);
		double exeTime = 0;
		double time_lastRun;
		double oneSliceTime;
		while (aks != nullptr) {
			for (int kI = 0; kI < BENCHMARK_KERNEL_SLICE_LOOP_NUM; ++kI) {
				exeTime += benchmarkKernelSlices(dp, aks, queues);
			}
			exeTime /= BENCHMARK_KERNEL_SLICE_LOOP_NUM;
			if (sliceFactor == 0) {
				oneSliceTime = exeTime;
				time_lastRun = exeTime;
			}
			// break if the time increases too much, 20% higher
//	  debug_printf(__func__, "Slice factor %d, one slice time: %lf, last run time:%lf, current execution time: %lf\n",
//				   sliceFactor, oneSliceTime, time_lastRun, exeTime);
			if ((exeTime - oneSliceTime) / oneSliceTime > 0.5
					|| (exeTime - time_lastRun) / time_lastRun > 0.2) { break; }
			time_lastRun = exeTime;
			aks = edcl_->slicingKernelPow2_extra_equal(kernel, ++sliceFactor);
		}
		// while loop breaks if time increases too much or kernel can't be further sliced.
		// either case sliceFactor = sliceFactor - 1;
		dp.NC = --sliceFactor;
	}

	void findNCKernelBased(EDKernel kernel, DP &dp) {
//		debug_printf(__func__, "\n");
		int sliceFactor = 0;
		KernelExecutionThread *execution_thread;
		bool cannotBeSlices = false;
		double time_lastRun;
		double oneSliceTime;
		while (!cannotBeSlices) {
			double avgTime = 0;
			dp.NC = sliceFactor;
			execution_thread = new KernelExecutionThread(edcl_, rm, dp, kernel);
			for (int i = 0; i < BENCHMARK_KERNEL_SLICE_LOOP_NUM; ++i) {
//				debug_printf(__func__, "[%d] Launching thread %s...\n", i, execution_thread->threadDescription.c_str());
				execution_thread->initThread();
				if (startKernelExecutionThread(execution_thread)) {
//					debug_printf(__func__, "Execution started, sliceFactor: %d, num chunks: %d. DP: ",
//											 sliceFactor, execution_thread->chunks->kernels.size());
//					execution_thread->dp.printDP();
					joinKernelExecutionThread(execution_thread);
					avgTime += execution_thread->activeTime;
//					execution_thread->printActiveTime();
//					debug_printf(__func__, "[%d] Thread %s joined! avgTime: %f\n",
//											 i, execution_thread->threadDescription.c_str(), avgTime);
				} else {
					cannotBeSlices = true; // can't be further sliced
//					debug_printf(__func__, "can't be further sliced\n");
					break;
				}
			}
			delete execution_thread;
			if (cannotBeSlices) break;
			avgTime /= BENCHMARK_KERNEL_SLICE_LOOP_NUM;
			if (sliceFactor == 0) {
				oneSliceTime = avgTime;
				time_lastRun = avgTime;
			}
			if ((avgTime - oneSliceTime) / oneSliceTime > 0.1 || (avgTime - time_lastRun) / time_lastRun > 0.1) { break; }
//			debug_printf(__func__, "sliceFactor: %d, oneSliceTime: %f, avgTime: %f, time_lastRun: %f\n",
//									 sliceFactor, oneSliceTime, avgTime, time_lastRun);
			sliceFactor++;
			time_lastRun = avgTime;
		}
		dp.NC = --sliceFactor;
		assert(dp.NC >= 0);
	}

	void setDPKernelBased(EDKernel kernel, DP &dp) {
//		debug_printf(__func__, "\n");
		float avgET = 0;
		KernelExecutionThread thread(edcl_, rm, dp, kernel);
		// Test run:
		if (!startKernelExecutionThread(&thread)) {
			err_printf(__func__, __LINE__, "%s starts failed! DP: \n");
			thread.dp.printDP();
			exit(EXIT_FAILURE);
		}
		joinKernelExecutionThread(&thread);
		// if no little core is used, set little core to lowest freq
		if (dp.nl < 0) rm->setCPUFrequency(0, rm->getHeteroCPUFrequencies()[0].size() - 1);
//		debug_printf(__func__, "Current L freq: %d\n", rm->getCPUCurrentFrequency(0));
		int numLoop = MAX(BENCHMARK_KERNEL_SLICE_LOOP_NUM, ceil(SET_DP_TIME / thread.activeTime));
		std::vector<float> current;
		current.reserve(100000);
		for (int i = 0; i < numLoop; ++i) {
			thread.initThread();
			rm->startCurrentRecording(&current);
			if (!startKernelExecutionThread(&thread)) {
				err_printf(__func__, __LINE__, "%s starts failed! DP: \n");
				thread.dp.printDP();
				exit(EXIT_FAILURE);
			}
			joinKernelExecutionThread(&thread);
			rm->stopCurrentRecording();
			avgET += (float)thread.activeTime;
		}
		float avgI = std::accumulate(current.begin(), current.end(), 0.0f) / current.size();
		avgET /= (float)numLoop;
		assert(avgI > 0);
		calculateEC(avgI, avgET, dp);
		dp.Prf = 1.0f / avgET;
		dp.CK = thread.dp.CK;
	}

	void genDPKernelBased(EDKernel kernel, DP &dp, std::vector<DP> &DPs) {
//		debug_printf(__func__, "\n");
		std::vector<float> prfs{0, 0, 0};
		getPerf(dp, prfs, DPs);
		float sumPrfs = 0;
		for (int kI = 0; kI < 3; ++kI) { sumPrfs += prfs[kI]; }
		for (int i = 0; i < 3; ++i) { dp.CK_ratio[i] = prfs[i] * (1.0 / sumPrfs); }
		findNCKernelBased(kernel, dp);
		setDPKernelBased(kernel, dp);
	}

// generate a DP, device based execution
	__unused void genDP(EDKernel kernel, DP &dp, std::vector<DP> &DPs) {
		// set exe queues and find proportion(determine CK_ratio)
		// 0: Little, 1: Big, 2: GPU
		std::vector<float> prfs{0, 0, 0};
		getPerf(dp, prfs, DPs);
		float sumPrfs = 0;
		for (int kI = 0; kI < 3; ++kI) { sumPrfs += prfs[kI]; }
		for (int i = 0; i < 3; ++i) { dp.CK_ratio[i] = prfs[i] * (1.0 / sumPrfs); }
		std::vector<EDQueue> queues{nullptr, nullptr, nullptr};
//	if (dp.nl >= 0) queues[0] = (rm->createQueueLittleCore(dp.nl, edcl_));
//	if (dp.nb >= 0) queues[1] = (rm->createQueueBigCore(dp.nb, edcl_));
//	if (dp.ng >= 0) queues[2] = (edcl_->createDeviceCmdQueueProfilingEnabled(GPUQueue));
		queues[0] = CHOOSE(dp.nl >= 0, rm->createQueueLittleCore(dp.nl, edcl_), nullptr);
		queues[1] = CHOOSE(dp.nb >= 0, rm->createQueueBigCore(dp.nb, edcl_), nullptr);
		queues[2] = CHOOSE(dp.ng >= 0, edcl_->createDeviceCmdQueueProfilingEnabled(GPUQueue), nullptr);

		// set frequency and execute to find the right NC
		rm->setCPUFrequency(0, dp.fl);
		rm->setCPUFrequency(rm->getNumCPU() - 1, dp.fb);
		rm->setGPUFrequency(dp.fg);
		// determine total number of chunks NC
		findNC(kernel, dp, queues);

		auto aks = edcl_->slicingKernelPow2_extra_equal(kernel, dp.NC);
		setDP(dp, aks, queues);

		for (auto &q: queues) {
			if (q != nullptr)
				edcl_->releaseEDQueue(q);
		}
	}

	static DP *findDP(std::vector<DP> &DPs, int nL, int nB, int nG, int fL, int fB, int fG) {
		for (auto &record : DPs) {
			if (nL == record.nl && nB == record.nb && nG == record.ng
					&& fL == record.fl && fB == record.fb && fG == record.fg)
				return &record;
		}
		return nullptr;
	}

	void genDPs(EDKernel k, std::vector<DP> &DPs, const char *csvFilename) { // generate DPs for a kernel
		time_t t = time(nullptr);
		struct tm tm = *localtime(&t);
		std::string filename = "NEW_" +
				std::to_string(tm.tm_year + 1900) +
				std::to_string(tm.tm_mon + 1) +
				std::to_string(tm.tm_mday) +
				std::to_string(tm.tm_hour) + std::to_string(tm.tm_min) + std::to_string(tm.tm_sec) + csvFilename;
		FILE *newDPRecord = fopen(filename.c_str(), "w");
		FILE *csvFile = fopen(csvFilename, "a");

		auto freq_L = rm->getHeteroCPUFrequencies()[0];
		auto freq_B = rm->getHeteroCPUFrequencies()[1];
		auto freq_G = rm->getGPUFrequencies();
		int t_nfL = freq_L.size(); // total num little core frequency
		int t_nfB = freq_B.size();
		int t_nfG = freq_G.size();
		double timer;

		for (int nL = -1; nL < 3; ++nL) {
			for (int nB = -1; nB < 3; ++nB) {
				for (int nG = -1; nG < 1; nG++) {
					for (int nfL = 0; nfL < t_nfL; ++nfL) {
						for (int nfB = 0; nfB < t_nfB; ++nfB) {
							for (int nfG = 0; nfG < t_nfG; ++nfG) {
								// DEBUG loops
//	for (int nL = -1; nL < 3; ++nL) {
//	  for (int nB = -1; nB < 3; ++nB) {
//		for (int nG = -1; nG < 1; nG++) {
//		  for (int nfL = 0; nfL < 1; ++nfL) {
//			for (int nfB = 0; nfB < 1; ++nfB) {
//			  for (int nfG = 0; nfG < 1; ++nfG) {
								// END DEBUG loops
								if (nB == -1 && nL == -1 && nG == -1) { // no need to loop through frequencies but no device is used
									nfL = t_nfL;
									nfG = t_nfG;
									nfB = t_nfB;
									continue;
								}
								// Skip generating redundant dps.
								// If there's only one device, there's no need to loop through frequencies for other devices
								if (nB != -1 && nL == -1 && nG == -1) {
									nfL = t_nfL;
									nfG = t_nfG;
								}
								if (nL != -1 && nB == -1 && nG == -1) {
									nfB = t_nfB;
									nfG = t_nfG;
								}
								if (nG != -1 && nL == -1 && nB == -1) {
									nfL = t_nfL;
									nfB = t_nfB;
								}
								// If there're only 2 devices used, there's no need to loop through frequencies for the other device
								if (nB != -1 && nL != -1 && nG == -1) { nfG = t_nfG; }
								if (nB != -1 && nL == -1 && nG != -1) { nfL = t_nfL; }
								if (nB == -1 && nL != -1 && nG != -1) { nfB = t_nfB; }

								if (findDP(DPs, nL, nB, nG, nfL, nfB, nfG) == nullptr) {
									DP dp;
									dp.DP_ID = genKernelID(k);
									dp.nl = nL;
									dp.nb = nB;
									dp.ng = nG;
									dp.fl = nfL;
									dp.fb = nfB;
									dp.fg = nfG;
									// DEBUG nl: 0, nb: 2, ng: 0, fl: 0, fb: 1, fg: 3
//				if (dp.nl == -1
//					&& dp.nb == 0
//					&& dp.ng == 0
//					&& dp.fl == 0
//					&& dp.fb == 0
//					&& dp.fg == 0) {
//				  debug_printf(__func__, " Hit bug point\n");
//				}
									debug_printf("Generating DP", "%s: nl: %d, nb: %d, ng: %d, fl: %d, fb: %d, fg: %d\n",
															 dp.DP_ID.c_str(), dp.nl, dp.nb, dp.ng, dp.fl, dp.fb, dp.fg);
									timer = getCurrentTime();
									// END DEBUG
//									genDP(k, dp, DPs);
									genDPKernelBased(k, dp, DPs);
//				sleep(1);
									// Print to keep track of progress
									printDPVertically(dp, rm);
									debug_printf("One DP time: ", "%lf\n", getCurrentTime() - timer);
									DPs.push_back(dp);
									writeDP(dp, newDPRecord);
									writeDP(dp, csvFile);
								}
							}
						}
					}
				}
			}
		}
	}

	static void printDPHead() {
		flushed_printf("DP_ID,");
		flushed_printf("Prf,");
		flushed_printf("EC,");
		flushed_printf("NC,");
		flushed_printf("CK_ratio,");
		flushed_printf("CK,");
		flushed_printf("nl,");
		flushed_printf("fl,");
		flushed_printf("nb,");
		flushed_printf("fb,");
		flushed_printf("ng,");
		flushed_printf("fg\n");
	}

	static void printDP(DP &dp) {
		flushed_printf("%s,", dp.DP_ID.c_str());
		flushed_printf("%f,", dp.Prf);
		flushed_printf("%f,", dp.EC);
		flushed_printf("%d,", dp.NC);

		flushed_printf("[ ");
		for (auto ratio:dp.CK_ratio) flushed_printf("%f ", ratio);
		flushed_printf("],");

		flushed_printf("[ ");
		for (auto ck:dp.CK) flushed_printf("%d ", ck);
		flushed_printf("],");

		flushed_printf("%d,", dp.nl);
		flushed_printf("%d,", dp.fl);
		flushed_printf("%d,", dp.nb);
		flushed_printf("%d,", dp.fb);
		flushed_printf("%d,", dp.ng);
		flushed_printf("%d\n", dp.fg);
	}

	static void printDPVertically(DP &dp, ResourceManager RM) {
		flushed_printf("DP_ID: %s\n", dp.DP_ID.c_str());
		flushed_printf("Prf: %f\n", dp.Prf);
		flushed_printf("EC: %f\n", dp.EC);
		flushed_printf("NC: %d\n", dp.NC);
		flushed_printf("CK: ");
		for (int kI = 0; kI < dp.CK.size(); ++kI) {
			flushed_printf("(%.3f, %d) ", dp.CK_ratio[kI], dp.CK[kI]);
		}
		flushed_printf("\n");
		if (dp.nl >= 0) {
			flushed_printf("nl: %d\n", dp.nl);
			flushed_printf("fl: %d\n", RM->getHeteroCPUFrequencies()[0].at(dp.fl));
		}
		if (dp.nb >= 0) {
			flushed_printf("nb: %d\n", dp.nb);
			flushed_printf("fb: %d\n", RM->getHeteroCPUFrequencies()[1].at(dp.fb));
		}
		if (dp.ng >= 0) {
			flushed_printf("ng: %d\n", dp.ng);
			flushed_printf("fg: %d\n", RM->getGPUFrequencies().at(dp.fg));
		}
	}

	static void writeDPHead(FILE *csvFile) {
		fprintf(csvFile, "DP_ID,");
		fprintf(csvFile, "Prf,");
		fprintf(csvFile, "EC,");
		fprintf(csvFile, "NC,");
		fprintf(csvFile, "CK_ratio,");
		fprintf(csvFile, "CK,");
		fprintf(csvFile, "nl,");
		fprintf(csvFile, "fl,");
		fprintf(csvFile, "nb,");
		fprintf(csvFile, "fb,");
		fprintf(csvFile, "ng,");
		fprintf(csvFile, "fg\n");
		fflush(csvFile);
	}

	static void writeDP(DP &dp, FILE *csvFile) {
		fprintf(csvFile, "%s,", dp.DP_ID.c_str());
		fprintf(csvFile, "%f,", dp.Prf);
		fprintf(csvFile, "%f,", dp.EC);
		fprintf(csvFile, "%d,", dp.NC);
		fprintf(csvFile, "[ ");
		for (auto ratio:dp.CK_ratio) fprintf(csvFile, "%f ", ratio);
		fprintf(csvFile, "],");

		fprintf(csvFile, "[ ");
		for (auto ck:dp.CK) fprintf(csvFile, "%lu ", ck);
		fprintf(csvFile, "],");

		fprintf(csvFile, "%d,", dp.nl);
		fprintf(csvFile, "%d,", dp.fl);
		fprintf(csvFile, "%d,", dp.nb);
		fprintf(csvFile, "%d,", dp.fb);
		fprintf(csvFile, "%d,", dp.ng);
		fprintf(csvFile, "%d\n", dp.fg);
		fflush(csvFile);
	}

	static void printTDP(TDP_MAP &TDP_in) {
		CA2020::printDPHead();
		for (auto &dp_p : TDP_in) {
			auto dp_v = dp_p.second;
			for (auto &dp : dp_v) {
				CA2020::printDP(dp);
			}
		}
	}

	static void writeTDPToCSVFile(const char *filename, TDP_MAP &TDP_in) {
		FILE *csvFile = fopen(filename, "w");
		CA2020::writeDPHead(csvFile);
		for (auto &dp_p : TDP_in) {
			auto dp_v = dp_p.second;
			for (auto &dp : dp_v) {
				CA2020::writeDP(dp, csvFile);
			}
		}
		fclose(csvFile);
	}

	static int readTDPFromCSVFile(const char *filename, TDP_MAP &TDP_in) {
		try {
			io::CSVReader<12> csv_reader(filename);
			csv_reader.read_header(io::ignore_extra_column,
														 "DP_ID", "Prf", "EC", "NC", "CK_ratio", "CK", "nl", "fl", "nb", "fb", "ng", "fg");
			DP tmp;
			std::vector<DP> allDPs;
			std::string ck_ratio;
			std::string ck;
			while (csv_reader.read_row(tmp.DP_ID,
																 tmp.Prf,
																 tmp.EC,
																 tmp.NC,
																 ck_ratio,
																 ck,
																 tmp.nl,
																 tmp.fl,
																 tmp.nb,
																 tmp.fb,
																 tmp.ng,
																 tmp.fg)) {
				DP dp;
				dp.DP_ID = tmp.DP_ID;
				dp.Prf = tmp.Prf;
				dp.EC = tmp.EC;
				dp.NC = tmp.NC;
				std::vector<float> tmpRatio;
				string2Floats(ck_ratio, tmpRatio);
				assert(tmpRatio.size() == 3);
				dp.CK_ratio = tmpRatio;
				std::vector<size_t> tmpCK;
				string2Size_t(ck, tmpCK);
				assert(tmpCK.size() == 3);
				dp.CK = tmpCK;
				dp.nl = tmp.nl;
				dp.fl = tmp.fl;
				dp.nb = tmp.nb;
				dp.fb = tmp.fb;
				dp.ng = tmp.ng;
				dp.fg = tmp.fg;
				allDPs.push_back(dp);
			}

			std::vector<DP> *dp_vec = nullptr;
			for (auto &dp: allDPs) {
				if (TDP_in.count(dp.DP_ID) < 1) {
					delete dp_vec;
					dp_vec = new std::vector<DP>;
					TDP_in[dp.DP_ID] = *dp_vec;
				}
				TDP_in[dp.DP_ID].push_back(dp);
			}
			delete dp_vec;
		} catch (io::error::can_not_open_file &) {
			return -1;
		}
		return 0;
	}

	static TDP_MAP distillTDP(const char *OutFileName, TDP_MAP &TDP_in) {
		TDP_MAP distilled;
		for (auto &dp_p : TDP_in) {
			auto dp_v = dp_p.second;
			std::vector<int> removeDPIdx;
			for (int i = 0; i < dp_v.size(); ++i) {
				if (isInVector(removeDPIdx, i)) continue;
				auto DPi = dp_v[i];
				// keep one device used
				if ((DPi.nl >= 0 && DPi.nb < 0 && DPi.ng < 0) ||
						(DPi.nl < 0 && DPi.nb >= 0 && DPi.ng < 0) ||
						(DPi.nl < 0 && DPi.nb < 0 && DPi.ng >= 0))
					continue;
				for (int j = i + 1; j < dp_v.size(); ++j) {
					if (isInVector(removeDPIdx, j)) continue;
					auto DPj = dp_v[j];
					// keep one device used
					if ((DPj.nl >= 0 && DPj.nb < 0 && DPj.ng < 0) ||
							(DPj.nl < 0 && DPj.nb >= 0 && DPj.ng < 0) ||
							(DPj.nl < 0 && DPj.nb < 0 && DPj.ng >= 0))
						continue;

					/*
					 * if execution time and energy consumption of
					 * a point using higher number of cores are the same or smaller
					 * than that of a point using lower number of cores,
					 * then the former point is discarded. All
					 * */
					// if both b and l are used
					if (DPi.nl >= 0 && DPj.nl >= 0 && DPi.nb >= 0 && DPj.nb >= 0 && DPi.ng == DPj.ng) {
						// erase i
						if (!isInVector(removeDPIdx, i) &&
								(DPi.nl >= DPj.nl && DPi.nb >= DPj.nb) && DPi.Prf <= DPj.Prf && DPi.EC >= DPj.EC) {
//			  debug_printf("Remove i", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(i);
						}
						// erase j
						if (!isInVector(removeDPIdx, j) &&
								(DPj.nl >= DPi.nl && DPj.nb >= DPi.nb) && DPj.Prf <= DPi.Prf && DPj.EC >= DPi.EC) {
//			  debug_printf("Remove j", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(j);
						}
					} // if only L
					else if (DPi.nb < 0 && DPj.nb < 0 && DPi.nl >= 0 && DPj.nl >= 0 && DPi.ng == DPj.ng) {
						// erase i
						if (!isInVector(removeDPIdx, i) &&
								DPi.nl >= DPj.nl && DPi.Prf <= DPj.Prf && DPi.EC >= DPj.EC) {
//			  debug_printf("Remove i", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(i);
						}
						// erase j
						if (!isInVector(removeDPIdx, j) &&
								DPj.nl >= DPi.nl && DPj.Prf <= DPi.Prf && DPj.EC >= DPi.EC) {
//			  debug_printf("Remove j", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(j);
						}
					} // if only B
					else if (DPi.nl < 0 && DPj.nl < 0 && DPi.nb >= 0 && DPj.nb >= 0 && DPi.ng == DPj.ng) {
						// erase i
						if (!isInVector(removeDPIdx, i) &&
								DPi.nb >= DPj.nb && DPi.Prf <= DPj.Prf && DPi.EC >= DPj.EC) {
//			  debug_printf("Remove i", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(i);
						}
						// erase j
						if (!isInVector(removeDPIdx, j) &&
								DPj.nb >= DPi.nb && DPj.Prf <= DPi.Prf && DPj.EC >= DPi.EC) {
//			  debug_printf("Remove j", "NL(%d, %d), NB(%d, %d), Prf(%f, %f), EC(%f, %f)\n",
//						   DPi.nl, DPj.nl, DPi.nb, DPj.nb, DPi.Prf, DPj.Prf, DPi.EC, DPj.EC);
							removeDPIdx.push_back(j);
						}
					}
				}
			}
//	  for (int r: removeDPIdx) {
//		printf("%d ", r);
//	  }
			std::vector<DP> dp_v_distill;
			for (int kI = 0; kI < dp_v.size(); ++kI) {
				bool distill = false;
				for (int r: removeDPIdx) {
					if (kI == r) distill = true;
				}
				if (!distill) dp_v_distill.push_back(dp_v[kI]);
			}
			std::sort(dp_v_distill.begin(), dp_v_distill.end(), std::greater<DP>());
			distilled[dp_p.first] = dp_v_distill;
			debug_printf(__func__, "%s: %d DPs removed from %d DPs, remain %d DPs.\n",
									 dp_p.first.c_str(), removeDPIdx.size(), dp_v.size(), dp_v_distill.size());
		}
		writeTDPToCSVFile(OutFileName, distilled);
		return distilled;
	}

	void printCDP(std::vector<DP> &CDP) const {
		int Nl = 0;
		int Nb = 0;
		float totalET = 0;
		float totalEC = 0;
		float ECETProduct = 0;
		for (auto &dp:CDP) {
			float fl = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f, 0);
			float fb = CHOOSE(dp.nb >= 0, (float)(rm->getHeteroCPUFrequencies()[1].at(dp.fb)) / 1000000.0f, 0);
			float fg = CHOOSE(dp.ng >= 0, (float)(rm->getGPUFrequencies().at(dp.fg)) / 1000.0f, 0);
			int nl = CHOOSE(dp.nl >= 0, pow2(dp.nl), 0);
			int nb = CHOOSE(dp.nb >= 0, pow2(dp.nb), 0);
			int ng = CHOOSE(dp.ng >= 0, pow2(dp.ng), 0);
			flushed_printf("\t");
			printDP(dp);
			flushed_printf("\t[%s] numL: %d, numB: %d, numG: %d. fl: %f, fb: %f, fg: %f. \n\t",
										 dp.DP_ID.c_str(), nl, nb, ng, fl, fb, fg);
			Nl += nl;
			Nb += nb;
			totalEC += dp.EC;
			totalET += 1 / dp.Prf;
			ECETProduct += dp.EC / dp.Prf;
		}
		totalEC *= 1000;
		ECETProduct *= 1000;
		flushed_printf("\tCDP info - [Core Count] L: %d, B: %d, EC: %f, ET: %f, Product: %f\n",
									 Nl, Nb, totalEC, totalET, ECETProduct);
	}

	friend class CA2020_UnitTest;
};

#endif //EDCL_INCLUDE_STG_CA2020_H_
