#ifndef EDCL_INCLUDE_EDCL_STRATEGY_H_
#define EDCL_INCLUDE_EDCL_STRATEGY_H_
#include "EDCL_ResourceManager.h"
#include "concurrentqueue.h"

class EDCL_Strategy {
 public: // variable
	double planingTime;
	double executionTime;
	explicit EDCL_Strategy(EDCL *edcl) {
		this->edcl_ = edcl;
		planingTime = 0;
		executionTime = 0;
	}

 protected: // variable
	EDCL *edcl_;

 protected: // functions don't expose to users but scheduler.
	std::vector<AtomicKernelSet> priorityStages{};
	// Offline scheduling, scheduling preparation
	virtual void plan(std::vector<AtomicKernelSet> &AKSs, ResourceManager ResourceManager) = 0;
	// Online scheduling
	virtual void execute(ResourceManager ResourceManager) = 0;
	// Releasing resources, free up tmp memory
	virtual void finish(ResourceManager ResourceManager) = 0;

	// add EDCL_Scheduler as a friend class for each strategy
	friend class EDCL_Scheduler;
};

inline void executeAtomicKernelSet(EDCL *edcl, AtomicKernelSet &kernelSet, EDQueue &queue) {
	double timer = getCurrentTime();
	for (auto &kernel : kernelSet->kernels) {
//	edcl->executeKernel(queue, kernel, 0, nullptr, kernel->kernelEvent->getCLEvent());
		edcl->executeKernel(queue,
												kernel,
												kernel->num_events_in_wait_list,
												kernel->event_wait_list,
												&kernel->kernelEvent);
	}
	for (auto &kernel : kernelSet->kernels) {
//    flushed_printf("Waiting for kernel %s...\n", kernel->name.c_str());
		kernel->waitForKernelEvent(queue);
	}
	timer = getCurrentTime() - timer;
	kernelSet->atomicExeTime = timer;
}
#endif //EDCL_INCLUDE_EDCL_STRATEGY_H_
