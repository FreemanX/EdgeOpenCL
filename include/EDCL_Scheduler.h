#ifndef EDCL_INCLUDE_EDCL_SCHEDULER_H_
#define EDCL_INCLUDE_EDCL_SCHEDULER_H_
#include "EDCL_Strategy.h"
#include "EDCL.h"
#include "Stg_Sequential.h"
#include <vector>
#include <mutex>

// Singleton
class EDCL_Scheduler {
 public: // variables
	std::vector<AtomicKernelSet> atomicKernelSets;

	// end public variables
 private: // variables
	EDCL *edcl_;
	ResourceManager resource_manager_;

// end private variables
 public: // functions
	static EDCL_Scheduler *getInstance(EDCL *edcl);
	EDCL_Scheduler(EDCL_Scheduler const &) = delete;
	void operator=(EDCL_Scheduler const &) = delete;
	ResourceManager getResourceManager() { return this->resource_manager_; }

	void submitKernels(EDKernel Kernel) {
		atomicKernelSets.push_back(edcl_->wrapEDKernelInAtomicKernelSet(Kernel));
	}

	template<typename ... EDKernels>
	void submitKernels(EDKernel Kernel, EDKernels ...Kernels) {
		submitKernels(Kernel);
		submitKernels(Kernels...);
	}

	void submitAtomicKernelSets(const AtomicKernelSet &kernelSet) {
		atomicKernelSets.push_back(kernelSet);
	}

	template<typename ... AtomicKernelSets>
	void submitAtomicKernelSets(const AtomicKernelSet &KernelSet, AtomicKernelSets ...KernelSets) {
		submitAtomicKernelSets(KernelSet);
		submitAtomicKernelSets(KernelSets...);
	}

	double executeKernels(EDCL_Strategy &strategy); // returns overall strategy execution time
	void clearKernels();
//  EDProgram createProgram(const char**ProgramSource); // reserved for kernel source transformation

	size_t getNumEDKernels() {
		size_t numKernels = 0;
		for (auto &ks : atomicKernelSets) {
			numKernels += ks->kernels.size();
		}
		return numKernels;
	}
	// end public functions
 private: // functions
	explicit EDCL_Scheduler(EDCL *edcl);

	int InitSchedulerEnv(std::vector<AtomicKernelSet> &priorityStage) {
		int numK = 0;
		std::vector<EDKernel> sortedKernels;
		// unmap
		for (auto &aks : atomicKernelSets) {
			for (auto &k : aks->kernels) {
				edcl_->unmapKernelBuffer(k);
				sortedKernels.push_back(k);
				numK++;
			}
		}

		// sort kernel based on priority in descending order
		std::sort(sortedKernels.begin(), sortedKernels.end(),
							[](const EDKernel &k1, const EDKernel &k2) -> bool {
								return k1->schedulePriority > k2->schedulePriority;
							});
		int lastPriority = 1;
		for (auto &k : sortedKernels) {
			if (lastPriority != k->schedulePriority) {
				priorityStage.push_back(edcl_->createAtomicKernelSet(1));
				lastPriority = k->schedulePriority;
			}
			priorityStage.back()->addKernel(k);
		}
		return numK;
	}

	~EDCL_Scheduler();

	friend class CA2020_UnitTest;
	friend class EdgeSchedule_UnitTest;

	// end private functions
};

// static function macros
#define GET_SCHEDULER(edcl) EDCL_Scheduler::getInstance(edcl)

#endif //EDCL_INCLUDE_EDCL_SCHEDULER_H_