#include <tar.h>
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#ifndef EDCL_INCLUDE_EDCL_H_
#define EDCL_INCLUDE_EDCL_H_
#include "custom_cl.h"
#include "CLRT.h"
#include <vector>
#include "BufferDef.h"
#include <memory>
#include <unordered_map>

enum QueueType { CPUQueue, GPUQueue, SDQueue };
typedef struct InternalEDCLBuffer_ *EDBuffer;
typedef struct InternalEDProgram_ *EDProgram;
typedef struct InternalEDKernel_ *EDKernel;
typedef struct InternalEDQueue_ *EDQueue;
typedef struct InternalEDUserEvent_ *EDUserEvent;
typedef struct InternalUserEventList_ *EDUserEventList;

typedef struct InternalEDProgram_ : public Trackable {
	CLRT_Program cpu_program{};
	CLRT_Program gpu_program{};
	__unused const char **sourceCode{};
} InternalEDProgram_;

class InternalEDQueue_ : public Trackable {
 public:
	CLRT_Queue clrt_queue{};
	cl_command_queue queue{};
	QueueType queueType{};
	// CPU affinity record
	u_int16_t cpuSet = 0;
	uint numCU{};
	uint collectionIdx{};
	uint coreStart{};
	CLRT *clrt; // which opencl runtime it belongs to
 public:
	explicit InternalEDQueue_(CLRT *CLRT) {
		clrt = CLRT;
	}
	void finish() const {
		CLRT_ERR(clFinish(queue, clrt->oclLibHandle_));
	}
	void *getSoHandle() const { return clrt->oclLibHandle_; }
};

class InternalEDUserEvent_ {
 public:
	cl_event gpuEvent{};
	cl_event cpuEvent{};
	InternalEDUserEvent_(CLRT *cpuRT, CLRT *gpuRT) {
		cpuRT_ = cpuRT;
		gpuRT_ = gpuRT;
		initEvents();
	}

	void initEvents() {
		releaseEvent(); // if exist previous event
		this->gpuEvent = gpuRT_->createUserEvent();
		this->cpuEvent = cpuRT_->createUserEvent();
	}

	void setEventComplete() {
		cpuRT_->setUserEventComplete(cpuEvent);
		gpuRT_->setUserEventComplete(gpuEvent);
	}

	void releaseEvent() {
		gpuRT_->releaseEvent(gpuEvent);
		cpuRT_->releaseEvent(cpuEvent);
	}

	~InternalEDUserEvent_() {
		releaseEvent();
	}
 private:
	CLRT *cpuRT_;
	CLRT *gpuRT_;
};

typedef struct InternalUserEventList_ {
	std::vector<EDUserEvent> EDUserEvents;

	void addEvent(EDUserEvent user_event) {
		EDUserEvents.push_back(user_event);
	}

	void setComplete() {
		for (auto &event : EDUserEvents)
			event->setEventComplete();
	}

	__unused void initEvents() {
		for (auto &event : EDUserEvents)
			event->initEvents();
	}

	void createEventArray(QueueType queue_type, cl_event *eventList, uint &numEvents) {
		numEvents = EDUserEvents.size();
		eventList = static_cast<cl_event *>(calloc(numEvents, sizeof(cl_event)));
		for (int i = 0; i < numEvents; i++)
			queue_type == GPUQueue ? eventList[i] = EDUserEvents[i]->gpuEvent : eventList[i] = EDUserEvents[i]->cpuEvent;
	}

	~InternalUserEventList_() {
		for (auto &event : EDUserEvents)
			delete event;
	}

} InternalUserEventList_;

struct InternalAtomicKernelSet_ : public Trackable {
	double atomicExeTime;
	std::vector<EDKernel> kernels;

	explicit InternalAtomicKernelSet_(uint numKernels) {
		kernels.reserve(numKernels);
		atomicExeTime = 0;
	}

	uint getNumKernels() const {
		return kernels.size();
	}

	void addKernel(EDKernel kernel) {
		kernels.push_back(kernel);
	}

	template<typename ...EDKernels>
	void addKernel(EDKernel kernel, EDKernels ...edkernels) {
		addKernel(kernel);
		addKernel(edkernels...);
	}
};
typedef std::shared_ptr<InternalAtomicKernelSet_> AtomicKernelSet;

// Edge Device openCL Runtime
class EDCL {
 public:
	CLRT *cpuRT_;
	CLRT *gpuRT_;

 private:
	std::unordered_map<int, EDBuffer> buffers_;
	std::unordered_map<int, EDProgram> programs_;
	std::unordered_map<int, EDQueue> queues_;
	std::unordered_map<int, EDKernel> kernels_;
	std::unordered_map<int, uint> kernelSliceTracker; // kernelID, num slices

 public:
	EDCL();

	EDBuffer createBuffer(size_t nByte);

	EDProgram createProgram(const char **ProgramSource);

	EDProgram createProgramWithOptions(const char **ProgramSource, const char *Options);

	EDKernel createKernel(EDProgram Program, std::string Name, uint NumArgs);

	EDKernel createKernelCopy(EDKernel kernel);

	AtomicKernelSet createAtomicKernelSet(uint numKernels);

	AtomicKernelSet wrapEDKernelInAtomicKernelSet(EDKernel Kernel);

	__unused AtomicKernelSet wrapEDKernelsInAtomicKernelSet(const std::vector<EDKernel> &Kernels);

	AtomicKernelSet slicingKernelPow2(EDKernel k, uint sliceFactor);

	AtomicKernelSet slicingKernelPow2_extra(EDKernel k, uint sliceFactor);

	AtomicKernelSet slicingKernelPow2_extra_equal(EDKernel k, uint sliceFactor);

	// TODO: unify queue request interface
	/*
	 * Log2numCU: How many cores needed?
	 * 	assume queue can only occupy 2^Log2numCU cores.
	 * 	Log2numCU can be 0, 1, 2, 3 for 8-core SoC
	 *
	 * SDIndex: Sub-Device Index. The SDIndex th sub-device
	 *	Each sub-device is bind to specific CPU cores.
	 *	SDIndex indicates which cores to use for the queue
	 *
	 * SDIndex indicating which sub-device the queue is associated
	 * E.g. for 4B + 4L core SoC:
	 * 	if Log2num = 0, then numCU=1 and SDIndex can be [0-7]
	 * 		CPU Affinity: <SDIndex>: <cpu id range>
	 * 		0: cpu[0], 1: cpu[1], 2: cpu[2] ...
	 * 	if Log2num = 1, then numCU=2 and SDIndex can be [0-3]
	 * 		0: cpu[0,1], 1: cpu[2,3], 2: cpu[4,5] ...
	 * 	if Log2num = 2, then numCU=4 and SDIndex can be [0-1]
	 * 		0: cpu[0,3], 1: cpu[4,7]
	 * 	if Log2num = 3, then numCU=8 and SDIndex can be [0]
	 * 		0: cpu[0,8]
	 * */
	// TODO: should not be Log2 of num of CUs. should be 1, 2, 4, 6...; e.g. 6 Little cores
	EDQueue createSubDeviceExeQueue(uint Log2numCU, uint SDIndex);

	// TODO: unify queue request interface
	EDQueue createDeviceCmdQueue(QueueType Type, cl_command_queue_properties Properties);

	// TODO: unify queue request interface
	EDQueue createDeviceCmdQueueProfilingEnabled(QueueType Type) {
		return createDeviceCmdQueue(Type, CL_QUEUE_PROFILING_ENABLE);
	}

	void confirmExeEnv(EDQueue Queue, EDKernel Kernel) const;

	__unused  void confirmAKSEnv(EDQueue Queue, const AtomicKernelSet &aks) const;

	void executeKernel(EDQueue Queue, EDKernel Kernel,
										 unsigned int NumEventsWait,
										 const cl_event *EventWaitList,
										 cl_event *event) const;

	void executeKernelBufferBlocked(EDQueue Queue, EDKernel Kernel,
																	unsigned int NumEventsWait,
																	const cl_event *EventWaitList,
																	cl_event *event);

	void unmapKernelBuffer(EDKernel Kernel);

	void mapKernelBuffer(EDKernel Kernel);

	void confirmExeEnvAndExecuteSingleKernel(EDQueue Queue,
																					 EDKernel Kernel,
																					 unsigned int NumEventsWait,
																					 const cl_event *EventWaitList,
																					 cl_event *event) {
		confirmExeEnv(Queue, Kernel);
		executeKernel(Queue, Kernel, NumEventsWait, EventWaitList, event);
		mapKernelBuffer(Kernel);
	}

	EDUserEvent createUserEvent() const;

	void releaseEDBuffer(EDBuffer Buffer);

	void releaseEDProgram(EDProgram Program);

	void releaseEDQueue(EDQueue Queue);

	uint getNumKernelSlices(int kernelID) {
		uint numSlice = 0;
		if (kernelSliceTracker.count(kernelID) > 0) numSlice = kernelSliceTracker[kernelID];
		return numSlice;
	}

	void resetNumKernelSlices(int kernelID) {
		if (kernelSliceTracker.count(kernelID) > 0) kernelSliceTracker[kernelID] = 0;
	}

	~EDCL();

 private:
	void mapEDBuffer(EDBuffer Buffer) const;

	void unmapEDBuffer(EDBuffer Buffer) const;

}; // end class

static int kernelCnt = 0;
class InternalEDKernel_ {
 public:
	InternalEDKernel_(std::string Name, uint NumArgs, EDProgram Program) {
		numArgs = NumArgs;
		name = std::move(Name);
		argTypeList = static_cast<BufferType *>(calloc(NumArgs, sizeof(BufferType)));
		argObjSizeList = static_cast<size_t *>(calloc(NumArgs, sizeof(size_t)));
		argPtrList = static_cast<void **>(calloc(NumArgs, sizeof(void *)));
		waitEventList = new InternalUserEventList_;
		subscriberEventList_ = new InternalUserEventList_;
		program_ = Program;
		kernelID = kernelCnt++;
	}

	void releaseKernelResource() {
		numArgs = -1;
		index = -1;
		workDim = -1;
		delete argTypeList;
		delete argObjSizeList;
		delete argPtrList;
		argTypeList = nullptr;
		argObjSizeList = nullptr;
		argPtrList = nullptr;
		delete waitEventList;
		waitEventList = nullptr;
		delete event_wait_list;
		event_wait_list = nullptr;
		delete subscriberEventList_;
		subscriberEventList_ = nullptr;
		delete[] globalSize;
		delete[] globalSizeOrg;
		delete[] localSize;
		delete[] offset;
		globalSize = nullptr;
		localSize = nullptr;
		offset = nullptr;
	}

	~InternalEDKernel_() {
		releaseKernelResource();
	}

	void printInfo() const {
		flushed_printf("Kernel %s info: offset(%d, %d, %d), global(%d, %d, %d), local(%d, %d, %d). Priority: %d\n",
									 name.c_str(),
									 CHOOSE(offset != nullptr, offset[0], 0),
									 CHOOSE(offset != nullptr && workDim > 1, offset[1], 0),
									 CHOOSE(offset != nullptr && workDim > 2, offset[2], 0),
									 globalSize[0],
									 CHOOSE(workDim > 1, globalSize[1], 0),
									 CHOOSE(workDim > 2, globalSize[2], 0),
									 localSize[0],
									 CHOOSE(workDim > 1, localSize[1], 0),
									 CHOOSE(workDim > 2, localSize[2], 0),
									 schedulePriority);
	}

	int kernelID;
	EDProgram program_;
	std::string name;
	int numArgs;
	int index{};
	CLRT_Kernel cpu_kernel{};
	CLRT_Kernel gpu_kernel{};
	// args
	BufferType *argTypeList;
	size_t *argObjSizeList; // size of actual buffer, not memory object(sizeof(EDBuffer))
	void **argPtrList; // EDBuffer or Primitive type ptr
	// exe settings
	uint workDim{};
	cl_event kernelEvent{}; // event for profiling and make sure this kernel complete
	double kernelTime = 0;

	int schedulePriority = 0;
	std::vector<EDKernel> subscribers;
	EDUserEventList waitEventList;  // events this kernel should wait for
	EDUserEventList subscriberEventList_; // events this kernel will set COMPLETE
	uint num_events_in_wait_list = 0;
	cl_event *event_wait_list = nullptr;
	EDQueue executionQ{};

	void confirmExecutionQ(EDQueue Queue) { executionQ = Queue; }

	void confirmWaitList(QueueType queue_type) {
		waitEventList->createEventArray(queue_type, event_wait_list, num_events_in_wait_list);
	}

	void updatePriority(int SuggestPriority, std::vector<int> &path) {
//	debug_printf(__func__, "Updating %s, suggested: %d, current: %d\n", name.c_str(),
//				 SuggestPriority, schedulePriority);
		for (auto id : path) {
			if (id == this->kernelID) {
				err_printf("Find cycle!", -1, " Cycle find in DAG!\n");
				exit(-1);
			}
		}
		path.push_back(this->kernelID);
		if (SuggestPriority < this->schedulePriority) {
			schedulePriority = SuggestPriority;
			int newSuggest = schedulePriority - 1;
			for (auto k : subscribers) {
				k->updatePriority(newSuggest, path);
			}
		}
	}

	// Kernel A subscribe to Kernel B --> B.addSubscriber(A)
	// Kernel A will wait for the completion of execution of Kernel B
	void addSubscriber(EDCL *edcl, EDKernel Kernel) {
		EDUserEvent event = edcl->createUserEvent();
		this->subscriberEventList_->addEvent(event);
		Kernel->waitEventList->addEvent(event);
		subscribers.push_back(Kernel);
	}

	void waitForKernelEvent() {
		clWaitForEvents(1, &kernelEvent, executionQ->clrt->oclLibHandle_);
		// notify all its subscribers
		subscriberEventList_->setComplete();
		kernelTime = executionQ->clrt->getExeTime(&kernelEvent);
		clReleaseEvent(kernelEvent, executionQ->clrt->oclLibHandle_);
	}

	void waitForKernelEvent(EDQueue queue) {
		clWaitForEvents(1, &kernelEvent, queue->clrt->oclLibHandle_);
		// notify all its subscribers
		subscriberEventList_->setComplete();
		kernelTime = queue->clrt->getExeTime(&kernelEvent);
		clReleaseEvent(kernelEvent, queue->clrt->oclLibHandle_);
	}

	double getKernelEventTime() const {
		return kernelTime;
	}

	template<typename ... Args>
	void configKernel(uint WorkDim, size_t *GlobalSize, size_t *LocalSize, Args ...argPtrs) {
		setExeConfig(WorkDim, GlobalSize, LocalSize);
		setArgs(argPtrs...);
	}

	template<typename T, typename ... Args>
	void setArgs(T arg, Args ...args) {
		setArgs(arg);
		setArgs(args...);
	}

	void setLocalSize(const size_t *LocalSize) {
		if (LocalSize == nullptr) return;
		if (localSize == nullptr) localSize = new size_t[workDim];
		for (int kI = 0; kI < workDim; ++kI) {
			localSize[kI] = LocalSize[kI];
		}
	}

	void setGlobalSize(const size_t *GlobalSize) {
		if (GlobalSize == nullptr) return;
		if (globalSize == nullptr) globalSize = new size_t[workDim];
		for (int kI = 0; kI < workDim; ++kI) {
			globalSize[kI] = GlobalSize[kI];
		}
	}

	void setOffset(const size_t *Offset) {
		if (Offset == nullptr) return;
		if (offset == nullptr) offset = new size_t[workDim];
		for (int kI = 0; kI < workDim; ++kI) {
			this->offset[kI] = Offset[kI];
//	  debug_printf(__func__, "Offset[%d] = %d, this->offset[%d] = %d\n",
//				   kI, Offset[kI], kI, offset[kI]);
		}
//	debug_printf(__func__, "");
//	printInfo();
	}

	void setLocalSize(uint WorkDim, const size_t *LocalSize) {
		workDim = WorkDim;
		if (LocalSize == nullptr) return;
		if (localSize == nullptr) localSize = new size_t[WorkDim];
		for (int kI = 0; kI < WorkDim; ++kI) {
			localSize[kI] = LocalSize[kI];
		}
	}

	void setGlobalSize(uint WorkDim, const size_t *GlobalSize) {
		workDim = WorkDim;
		if (GlobalSize == nullptr) return;
		if (globalSize == nullptr) globalSize = new size_t[WorkDim];
		for (int kI = 0; kI < WorkDim; ++kI) {
			globalSize[kI] = GlobalSize[kI];
		}
	}

	__unused void setOffset(uint WorkDim, const size_t *Offset) {
		workDim = WorkDim;
		if (Offset == nullptr) return;
		if (offset == nullptr) offset = new size_t[WorkDim];
		for (int kI = 0; kI < WorkDim; ++kI) {
			offset[kI] = Offset[kI];
		}
	}

	void setGlobalSizeOrg(uint WorkDim, const size_t *GlobalSize) {
		workDim = WorkDim;
		if (GlobalSize == nullptr) return;
		if (globalSize == nullptr) globalSizeOrg = new size_t[WorkDim];
		for (int kI = 0; kI < WorkDim; ++kI) {
			globalSizeOrg[kI] = GlobalSize[kI];
		}
	}

	void setGlobalSizeOrg(const size_t *GlobalSize) {
		if (GlobalSize == nullptr) return;
		if (globalSize == nullptr) globalSizeOrg = new size_t[workDim];
		for (int kI = 0; kI < workDim; ++kI) {
			globalSizeOrg[kI] = GlobalSize[kI];
		}
	}

	size_t *getGlobalSize() { return this->globalSize; }
	size_t *getGlobalSizeOrg() { return this->globalSizeOrg; }
	size_t *getLocalSize() { return this->localSize; }
	size_t *getOffset() { return this->offset; }

	// calling this function will set the original GS
	void setExeConfig(uint WorkDim, size_t *GlobalSize, size_t *LocalSize) {
		if (offset == nullptr) offset = static_cast<size_t *>(calloc(WorkDim, sizeof(size_t)));
		if (globalSizeOrg == nullptr) globalSizeOrg = static_cast<size_t *>(calloc(WorkDim, sizeof(size_t)));
		if (globalSize == nullptr) globalSize = static_cast<size_t *>(calloc(WorkDim, sizeof(size_t)));
		if (LocalSize != nullptr && localSize == nullptr)
			localSize = static_cast<size_t *>(calloc(WorkDim, sizeof(size_t)));
		this->workDim = WorkDim;
		setGlobalSizeOrg(WorkDim, GlobalSize); // keep record of original kernel GS.
		setGlobalSize(WorkDim, GlobalSize);
		setLocalSize(WorkDim, LocalSize);
	}

	template<typename T>
	__unused void setArgAtIndex(T arg, uint idx) {
		size_t argSize = sizeof(*arg);
		if (argSize == sizeof(InternalEDCLBuffer_)) {
			auto edBuffer = (EDBuffer)arg;
			if (edBuffer->bufferType == ED_BUF) {
				this->argObjSizeList[idx] = edBuffer->sizeByte;
				this->argPtrList[idx] = arg;
				this->argTypeList[idx] = ED_BUF;
			}
		} else {
			this->argObjSizeList[idx] = argSize;
			this->argPtrList[idx] = arg;
			this->argTypeList[idx] = PRIMITIVE;
		}
	}

	bool operator<(const InternalEDKernel_ &kernel) const {
		return (schedulePriority < kernel.schedulePriority);
	}

	bool operator>(const InternalEDKernel_ &kernel) const {
		return (schedulePriority > kernel.schedulePriority);
	}

 private:
	size_t *offset = nullptr;
	size_t *globalSize = nullptr;
	size_t *globalSizeOrg = nullptr;
	size_t *localSize = nullptr;

	template<typename T>
	void setArgs(T arg) {
		uint idx = this->index % this->numArgs;
		size_t argSize = sizeof(*arg);
		if (argSize == sizeof(InternalEDCLBuffer_)) {
			auto edBuffer = (EDBuffer)arg;
			if (edBuffer->bufferType == ED_BUF) {
				this->argObjSizeList[idx] = edBuffer->sizeByte;
				this->argPtrList[idx] = arg;
				this->argTypeList[idx] = ED_BUF;
			}
		} else {
			this->argObjSizeList[idx] = argSize;
			this->argPtrList[idx] = arg;
			this->argTypeList[idx] = PRIMITIVE;
		}
		this->index++;
	}
};

#endif // EDCL_INCLUDE_EDCL_H_
