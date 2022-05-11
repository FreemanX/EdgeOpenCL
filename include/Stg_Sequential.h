#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#ifndef EDCL_INCLUDE_STG_SEQUENTIAL_H_
#define EDCL_INCLUDE_STG_SEQUENTIAL_H_
#include "EDCL_Strategy.h"

enum EXECUTE_DEVICE {
  B_SD1 = 0x1 << 0, // 1 big core sub-device, values are made-up
  B_SD2 = 0x2 << 1, // 2 big core sub-device
  B_SD4 = 0x3 << 2, // ...
  L_SD1 = 0x4 << 3, // 1 little core sub-device
  L_SD2 = 0x5 << 4, // ...
  L_SD4 = 0x5 << 5, // ...
  GPU = 0x6 << 6,
  CPU = 0x7 << 7,
};

inline char const *EXECUTE_DEVICE_ToString(EXECUTE_DEVICE device) {
  switch (device) {
	case CPU: {
	  return "CPU";
	}
	case GPU: {
	  return "GPU";
	}
	case B_SD1: {
	  return "1 CU Sub-device on Big core";
	}
	case B_SD2: {
	  return "2 CU Sub-device on Big core";
	}
	case B_SD4: {
	  return "4 CU Sub-device on Big core";
	}
	case L_SD1: {
	  return "1 CU Sub-device on Little core";
	}
	case L_SD2: {
	  return "2 CU Sub-device on Little core";
	}
	case L_SD4: {
	  return "4 CU Sub-device on Little core";
	}
	default: {
	  return "Unknown Device";
	}
  }
}

typedef struct ThreadWorkPack_ *ThreadWorkPack;
using namespace moodycamel;
struct ThreadWorkPack_ {
  EDCL *edcl_;
  EDQueue ed_queue_;
  ConcurrentQueue<AtomicKernelSet> *atomicKernelSetPool;
  ThreadWorkPack_(EDCL *edcl, EDQueue ed_queue, size_t NumKernels) {
	edcl_ = edcl;
	ed_queue_ = ed_queue;
	atomicKernelSetPool = new ConcurrentQueue<AtomicKernelSet>(NumKernels);
  }

  void enqueueAKS(const AtomicKernelSet& AKS) const {
	this->atomicKernelSetPool->enqueue(AKS);
  }

  ~ThreadWorkPack_() {
	delete atomicKernelSetPool;
  }
};

// PThread function
inline void *workerThread_Seq(void *WorkPack) {
  // Set kernel launching thread affinity
  setThreadAffinity(1);
  // Unpack work pack
  auto workPack = (ThreadWorkPack)WorkPack;
  auto edcl = workPack->edcl_;
  auto queue = workPack->ed_queue_;
  auto kernelSets = workPack->atomicKernelSetPool;
  // execute
  AtomicKernelSet kernelSet;
  while (kernelSets->try_dequeue(kernelSet)) {
	executeAtomicKernelSet(edcl, kernelSet, queue);
  }
  // finish
  queue->finish();
  pthread_exit(nullptr);
}

/* Strategy Implementations*/
class Sequential : public EDCL_Strategy {
 public:
  explicit Sequential(EDCL *edcl) : EDCL_Strategy(edcl) {
	execute_device_ = GPU;
  }

  void setExecutionDevice(EXECUTE_DEVICE Device) {
	this->execute_device_ = Device;
  }

 protected:
  EXECUTE_DEVICE execute_device_;
  EDQueue queue{};
  ThreadWorkPack workPack{};

  void plan(std::vector<AtomicKernelSet> &AKSs, ResourceManager resourceManager) override {
	planingTime = getCurrentTime();
	createCmdQueue();
	if (resourceManager->bindEDQueueCPU(queue) != 0) {
	  std::cerr << " Core binding failed!\n";
	  return;
	}

	workPack = new ThreadWorkPack_(edcl_, queue, AKSs.size());
	for (auto &kernelSet : AKSs) {
	  for (auto &kernel : kernelSet->kernels)
		edcl_->confirmExeEnv(queue, kernel);
	  workPack->atomicKernelSetPool->enqueue(kernelSet);
	}
	planingTime = getCurrentTime() - planingTime;
  }

  void execute(ResourceManager resourceManager) override {
	pthread_t pthread;
	executionTime = getCurrentTime();
	pthread_create(&pthread, nullptr, workerThread_Seq, workPack);
	pthread_join(pthread, nullptr);
	executionTime = getCurrentTime() - executionTime;
  }

  void finish(ResourceManager resourceManager) override {
	resourceManager->unbindEDQueueCPU(queue);
	delete workPack;
	edcl_->releaseEDQueue(queue);
  }

  void createCmdQueue() {
	switch (execute_device_) {
	  case CPU: {
		queue = edcl_->createDeviceCmdQueueProfilingEnabled(CPUQueue);
		break;
	  }
	  case GPU: {
		queue = edcl_->createDeviceCmdQueueProfilingEnabled(GPUQueue);
		break;
	  }
	  case B_SD1: {
		queue = edcl_->createSubDeviceExeQueue(0, 7);
		break;
	  }
	  case B_SD2: {
		queue = edcl_->createSubDeviceExeQueue(1, 3);
		break;
	  }
	  case B_SD4: {
		queue = edcl_->createSubDeviceExeQueue(2, 1);
		break;
	  }
	  case L_SD1: {
		queue = edcl_->createSubDeviceExeQueue(0, 3);
		break;
	  }
	  case L_SD2: {
		queue = edcl_->createSubDeviceExeQueue(1, 1);
		break;
	  }
	  case L_SD4: {
		queue = edcl_->createSubDeviceExeQueue(2, 0);
		break;
	  }
	  default: {
		std::cerr << "Unknown device, please reset\n";
		exit(-1);
	  }
	}

  }

  friend class EDCL_Scheduler;
};

class SequentialImproved : public Sequential {
 protected:
  void plan(std::vector<AtomicKernelSet> &AKSs, ResourceManager resourceManager) override {
	planingTime = getCurrentTime();
	createCmdQueue();
	if (resourceManager->bindEDQueueCPU(queue) != 0) {
	  std::cerr << " Core binding failed!\n";
	  return;
	}
	workPack = new ThreadWorkPack_(edcl_, queue, AKSs.size());
	AtomicKernelSet seqKernels = edcl_->createAtomicKernelSet(AKSs.size());
	for (auto &kernelSet : AKSs) {
	  for (auto &kernel : kernelSet->kernels)
		edcl_->confirmExeEnv(queue, kernel);
	  if (kernelSet->getNumKernels() < 2) {
		seqKernels->addKernel(kernelSet->kernels.back());
	  } else {
		workPack->enqueueAKS(seqKernels);
		seqKernels = edcl_->createAtomicKernelSet(AKSs.size() - seqKernels->getNumKernels());
		workPack->atomicKernelSetPool->enqueue(kernelSet);
	  }
	}
	if (seqKernels->getNumKernels() > 0) workPack->enqueueAKS(seqKernels);

	planingTime = getCurrentTime() - planingTime;
  }
 public:
  explicit SequentialImproved(EDCL *p_edcl) : Sequential(p_edcl) {}
  friend class EDCL_Scheduler;
};

#endif //EDCL_INCLUDE_STG_SEQUENTIAL_H_
