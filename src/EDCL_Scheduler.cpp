#include "EDCL_Scheduler.h"
#include "CPU_utils.h"

static std::mutex mutex_; // lock for thread-safe get instance
static EDCL_Scheduler *scheduler_;

EDCL_Scheduler *EDCL_Scheduler::getInstance(EDCL *edcl) {
	std::lock_guard<std::mutex> lock(mutex_);
	if (scheduler_ == nullptr) {
		scheduler_ = new EDCL_Scheduler(edcl);
	}
	return scheduler_;
}

EDCL_Scheduler::EDCL_Scheduler(EDCL *edcl) {
	this->edcl_ = edcl;
	this->resource_manager_ = new ResourceManager_();
}

EDCL_Scheduler::~EDCL_Scheduler() { clearKernels(); }

void EDCL_Scheduler::
clearKernels() {
	if (!atomicKernelSets.empty()) atomicKernelSets.clear();
}

double EDCL_Scheduler::
executeKernels(EDCL_Strategy &strategy) {
	if (atomicKernelSets.empty()) {
		err_printf(__func__, __LINE__, "No kernel to execute, please submit kernels\n");
		return -1;
	}

	std::vector<AtomicKernelSet> priorityStages;
	int numK = InitSchedulerEnv(priorityStages);

	strategy.plan(priorityStages, resource_manager_);

	strategy.execute(resource_manager_);

	strategy.finish(resource_manager_);

	debug_printf(__func__, "Num stages: %d, Total kernels: %d \n", priorityStages.size(), numK);

	// map
	for (auto &aks : atomicKernelSets) {
		for (auto &k : aks->kernels)
			edcl_->mapKernelBuffer(k);
	}

//  return strategy.executionTime + strategy.planingTime;
	return strategy.executionTime;
}
