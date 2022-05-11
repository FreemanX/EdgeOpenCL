#include <elf.h>
#pragma ide diagnostic ignored "readability-convert-member-functions-to-static"
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#pragma ide diagnostic ignored "OCDFAInspection"
#include "EDCL.h"
#include "utils.h"
#include <utility>
#include <CPU_utils.h>


// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()

EDCL::EDCL() {
	this->cpuRT_ = CLRT::createCPURunTime();
	this->gpuRT_ = CLRT::createGPURunTime();
}

EDBuffer EDCL::createBuffer(size_t nByte) {
	auto edBuffer = new InternalEDCLBuffer_(nByte);
	edBuffer->GPUzBuffer = gpuRT_->createZBuffer(nByte, nullptr);
	edBuffer->hostPtr = edBuffer->GPUzBuffer->hostPtr;
	edBuffer->CPUBuffer = cpuRT_->createBuffer(CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, nByte, edBuffer->hostPtr);
//  debug_printf(__func__, "%d\n", edBuffer->getTrackNum());
	this->buffers_[edBuffer->getTrackNum()] = edBuffer;
	return edBuffer;
}

EDProgram EDCL::createProgram(const char **ProgramSource) {
	auto edProgram = new InternalEDProgram_;
	edProgram->gpu_program = gpuRT_->createProgramFromSource(ProgramSource);
	edProgram->cpu_program = cpuRT_->createProgramFromSource(ProgramSource);
	this->programs_[edProgram->getTrackNum()] = edProgram;
	return edProgram;
}

EDProgram EDCL::createProgramWithOptions(const char **ProgramSource, const char *Options) {
	auto edProgram = new InternalEDProgram_;
	edProgram->cpu_program = cpuRT_->createProgramFromSourceWithOptions(ProgramSource, Options);
	edProgram->gpu_program = gpuRT_->createProgramFromSourceWithOptions(ProgramSource, Options);
	edProgram->sourceCode = ProgramSource;
	this->programs_[edProgram->getTrackNum()] = edProgram;
	return edProgram;
}

EDKernel EDCL::createKernel(EDProgram Program, std::string Name, uint NumArgs) {
	auto edKernel = new InternalEDKernel_(std::move(Name), NumArgs, Program);
	edKernel->cpu_kernel = cpuRT_->createKernel(Program->cpu_program, edKernel->name, edKernel->numArgs);
	edKernel->gpu_kernel = gpuRT_->createKernel(Program->gpu_program, edKernel->name, edKernel->numArgs);
	kernels_[edKernel->kernelID] = edKernel;
	return edKernel;
}

/*
 * Not copied but points to input kernel(new edkernel and kernel shared objects):
 * 	cpu_kernel, gpu_kernel
 *
 * Copied(new edkernel copies from kernel and keep its own copies):
 * 	name, numArgs, argTypeList, argObjSizeList, argPtrList, workDim, globalSize, localSize
 *
 * Ignored(init like creating new EDKernel):
 * 	index, offset, kernelEvent
 *
 * */
EDKernel EDCL::createKernelCopy(EDKernel kernel) {
	auto edKernel = new InternalEDKernel_(kernel->name, kernel->numArgs, kernel->program_);
//  edKernel->cpu_kernel = cpuRT_->createKernel(kernel->program_->cpu_program, edKernel->name, edKernel->numArgs);
//  edKernel->gpu_kernel = gpuRT_->createKernel(kernel->program_->gpu_program, edKernel->name, edKernel->numArgs);
	edKernel->cpu_kernel = kernel->cpu_kernel;
	edKernel->gpu_kernel = kernel->gpu_kernel;
	edKernel->kernelID = kernel->kernelID; // Copied kernel will have the same kernel ID
	if (kernelSliceTracker.count(kernel->kernelID) > 0) {
		kernelSliceTracker[kernel->kernelID] += 1;
	} else {
		kernelSliceTracker[kernel->kernelID] = 1;
	}

	for (int i = 0; i < kernel->numArgs; ++i) {
		edKernel->argTypeList[i] = kernel->argTypeList[i];
		edKernel->argObjSizeList[i] = kernel->argObjSizeList[i];
		edKernel->argPtrList[i] = kernel->argPtrList[i];
	}

	edKernel->setExeConfig(kernel->workDim, kernel->getGlobalSize(), kernel->getLocalSize());
	edKernel->setGlobalSizeOrg(kernel->getGlobalSizeOrg());

	return edKernel;
}

void EDCL::releaseEDBuffer(EDBuffer Buffer) {
	if (buffers_.count(Buffer->getTrackNum()) > 0) {
//	debug_printf(__func__, "%d\n", Buffer->getTrackNum());
		this->gpuRT_->releaseCLRT_zBuffer(Buffer->GPUzBuffer);
		this->cpuRT_->releaseCLRT_Buffer(Buffer->CPUBuffer);
		buffers_.erase(Buffer->getTrackNum());
		delete Buffer;
	}
}

void EDCL::releaseEDProgram(EDProgram Program) {
	if (programs_.count(Program->getTrackNum()) > 0) {
		this->cpuRT_->releaseCLRT_Program(Program->cpu_program);
		this->gpuRT_->releaseCLRT_Program(Program->gpu_program);
		programs_.erase(Program->getTrackNum());
		delete Program;
	}

}

EDCL::~EDCL() {
	std::vector<int> keys;

	for (auto &b: buffers_) { keys.push_back(b.first); }
	for (auto k :keys) { releaseEDBuffer(buffers_[k]); }
	keys.clear();

	for (auto &p: programs_) { keys.push_back(p.first); }
	for (auto &k: keys) { releaseEDProgram(programs_[k]); }
	keys.clear();

	for (auto &q: queues_) keys.push_back(q.first);
	for (auto &k: keys) releaseEDQueue(queues_[k]);

	delete cpuRT_;
	delete gpuRT_;
}

EDQueue EDCL::
createSubDeviceExeQueue(uint Log2numCU, uint SDIndex) {
	if (Log2numCU > cpuRT_->sub_deviceCollection.size()) {
		std::cerr << "NumCU " << exp2(Log2numCU) << " exceeds the total number of CPUs \n";
		exit(-1);
	}
	auto q = new InternalEDQueue_(cpuRT_);
	q->numCU = exp2(Log2numCU);
	q->cpuSet = 0;
	q->collectionIdx = Log2numCU;
	q->queueType = SDQueue;
	q->clrt_queue = cpuRT_->createSubDeviceExecutionProfilingQueue(Log2numCU, SDIndex);
	q->queue = q->clrt_queue->command_queue;
	q->coreStart = SDIndex * q->numCU;
	for (uint i = q->coreStart; i < q->coreStart + q->numCU; ++i)
		SET_BIT(q->cpuSet, i);
	queues_[q->getTrackNum()] = q;
	return q;
}

EDQueue EDCL::
createDeviceCmdQueue(QueueType Type, cl_command_queue_properties Properties) {
	if (Type == SDQueue) {
		std::cerr << "Use getSubDeviceExeQueue() to get SDQueue.\n";
		exit(-1);
	}
	EDQueue queue;
	CLRT *clrt;
	Type == CPUQueue ? clrt = cpuRT_ : clrt = gpuRT_;
	queue = new InternalEDQueue_(clrt);
	queue->queueType = Type;
	queue->clrt_queue = clrt->createCLExecutionQueue(Properties);
	queue->queue = queue->clrt_queue->command_queue;
	if (Type == CPUQueue) {
		queue->cpuSet = exp2(getNumCPUs()) - 1;
		queue->numCU = getNumCPUs();
		queue->coreStart = 0;
	} else if (Type == GPUQueue) {
		SET_BIT(queue->cpuSet, 1);
	}
	queues_[queue->getTrackNum()] = queue;
	return queue;
}

void EDCL::releaseEDQueue(EDQueue Queue) {
	if (queues_.count(Queue->getTrackNum()) > 0) {
		Queue->clrt->releaseCLRT_Queue(Queue->clrt_queue);
		queues_.erase(Queue->getTrackNum());
		delete Queue;
	}
}

void EDCL::
mapEDBuffer(EDBuffer Buffer) const {
	Buffer->hostPtr = gpuRT_->mapZBuffer(Buffer->GPUzBuffer);
}

void EDCL::
unmapEDBuffer(EDBuffer Buffer) const {
	gpuRT_->unmapZBuffer(Buffer->GPUzBuffer);
}

void EDCL::
unmapKernelBuffer(EDKernel Kernel) {
	for (int i = 0; i < Kernel->numArgs; ++i) {
		if (Kernel->argTypeList[i] == ED_BUF) {
			unmapEDBuffer((EDBuffer)Kernel->argPtrList[i]);
		}
	}
}

void EDCL::
mapKernelBuffer(EDKernel Kernel) {
	for (int i = 0; i < Kernel->numArgs; ++i) {
		if (Kernel->argTypeList[i] == ED_BUF) {
			mapEDBuffer((EDBuffer)Kernel->argPtrList[i]);
		}
	}
}

void EDCL::
confirmExeEnv(EDQueue Queue, EDKernel Kernel) const {
	if (Queue == nullptr) {
		err_printf(__func__, __LINE__, "Kernel %s is trying to bind with null queue!\n", Kernel->name.c_str());
	}
	if (Kernel == nullptr) {
		err_printf(__func__, __LINE__, "Trying to bind a Queue with null Kernel!\n");
	}
	Kernel->confirmWaitList(Queue->queueType);
	for (int i = 0; i < Kernel->numArgs; ++i) {
		switch (Kernel->argTypeList[i]) {
			case ED_BUF: {
				auto buffer = (EDBuffer)Kernel->argPtrList[i];
				if (Queue->queueType == GPUQueue) {
					gpuRT_->setKernelArgAtIndex(
							Kernel->gpu_kernel, i,
							buffer->GPUzBuffer->SIGNATURE,
							Kernel->argObjSizeList[i],
							buffer->GPUzBuffer);
				} else {
					cpuRT_->setKernelArgAtIndex(
							Kernel->cpu_kernel, i,
							buffer->CPUBuffer->SIGNATURE,
							Kernel->argObjSizeList[i],
							buffer->CPUBuffer);
				}
				// Once confirmed, EDBuffer won't be accessible via host ptr
//				unmapEDBuffer(buffer);
				break;
			}
			case PRIMITIVE: {
				CLRT *clrt_;
				Queue->queueType == GPUQueue ?
						clrt_ = gpuRT_ : clrt_ = cpuRT_;
				CLRT_Kernel kernel_;
				Queue->queueType == GPUQueue ?
						kernel_ = Kernel->gpu_kernel : kernel_ = Kernel->cpu_kernel;
				clrt_->setKernelArgAtIndex(
						kernel_,
						i,
						PRIMITIVE,
						Kernel->argObjSizeList[i],
						Kernel->argPtrList[i]);
				break;
			}
//	  case RELEASED: {
//		flushed_printf("Arg %d, Error!, released Buffer, size: %d\n", i, Kernel->argObjSizeList[i]);
//		exit(-1);
//	  }
			default: {
				flushed_printf("Kernel %s, ARG index %d: ", Kernel->name.c_str(), i);
				if (Kernel->argObjSizeList[i] < 1)
					flushed_printf("not set! size: %d\n", i, Kernel->argObjSizeList[i]);
				else
					flushed_printf("Error! Unknown Type! size: %d\n", i, Kernel->argObjSizeList[i]);
				exit(-1);
			}
		}
	}
}

void EDCL::
executeKernelBufferBlocked(EDQueue Queue,
													 EDKernel Kernel,
													 unsigned int NumEventsWait,
													 const cl_event *EventWaitList,
													 cl_event *event) {
	unmapKernelBuffer(Kernel);
	CLRT *clrt_;
	Queue->queueType == GPUQueue ? clrt_ = gpuRT_ : clrt_ = cpuRT_;
	CLRT_Kernel kernel_;
	Queue->queueType == GPUQueue ? kernel_ = Kernel->gpu_kernel : kernel_ = Kernel->cpu_kernel;
	clrt_->execKernel(
			Queue->queue, kernel_, Kernel->workDim, Kernel->offset,
			Kernel->globalSize, Kernel->localSize,
			NumEventsWait, EventWaitList, event);
	mapKernelBuffer(Kernel);
}

void EDCL::
executeKernel(EDQueue Queue,
							EDKernel Kernel,
							unsigned int NumEventsWait,
							const cl_event *EventWaitList,
							cl_event *event) const {
	CLRT *clrt_ = Queue->clrt;
	CLRT_Kernel kernel_;
	Queue->queueType == GPUQueue ? kernel_ = Kernel->gpu_kernel : kernel_ = Kernel->cpu_kernel;
	Kernel->confirmExecutionQ(Queue);
	clrt_->execKernel(Queue->queue, kernel_, Kernel->workDim, Kernel->offset,
										Kernel->globalSize, Kernel->localSize,
										NumEventsWait, EventWaitList, event);
}

AtomicKernelSet EDCL::
createAtomicKernelSet(uint numKernels) {
	auto atomic_kernel_set = std::make_shared<InternalAtomicKernelSet_>(numKernels);
	return atomic_kernel_set;
}

AtomicKernelSet EDCL::
wrapEDKernelInAtomicKernelSet(EDKernel Kernel) {
	AtomicKernelSet aks = createAtomicKernelSet(1);
	aks->addKernel(Kernel);
	return aks;
}

__unused AtomicKernelSet EDCL::
wrapEDKernelsInAtomicKernelSet(const std::vector<EDKernel> &Kernels) {
	AtomicKernelSet aks = createAtomicKernelSet(Kernels.size());
	for (auto &k : Kernels)
		aks->addKernel(k);
	return aks;
}

EDUserEvent EDCL::
createUserEvent() const {
	auto event = new InternalEDUserEvent_(cpuRT_, gpuRT_);
	return event;
}

AtomicKernelSet EDCL::
slicingKernelPow2(EDKernel k, uint sliceFactor) {
	if (k->localSize == nullptr) {
		PRINT_RED
		flushed_printf("Kernel %s local work size not set!\n", k->name.data());
		PRINT_DEFAULT
		return nullptr;
	}
	uint kDim = k->workDim; // kernel NDRange dimensions
	for (int kI = 0; kI < kDim; ++kI) {
		if (k->globalSize[kI] % k->localSize[kI] != 0) {
			PRINT_RED
			flushed_printf("Kernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n",
										 k->name.data(),
										 kI, k->globalSize[kI], k->localSize[kI]);
			PRINT_DEFAULT
			return nullptr;
		}
	}
	uint requiredNumSlices = pow2(sliceFactor);
	std::unique_ptr<size_t[]> numWGs = std::make_unique<size_t[]>(kDim);
	uint maxNumSlices = 1;
	for (int kI = 0; kI < kDim; ++kI) {
		numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
		maxNumSlices *= numWGs[kI];
	}
	if (maxNumSlices < requiredNumSlices) {
		PRINT_RED
		flushed_printf("Slice factor %d exceeds the kernel global size. ", sliceFactor);
		flushed_printf("Number of WGs({");
		for (int kI = 0; kI < kDim; ++kI)
			flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
		flushed_printf("}):%d. Required number of slices: %d.\n", maxNumSlices, requiredNumSlices);
		PRINT_DEFAULT
		return nullptr;
	}

	UNIQUE_PTR_SIZE_T(numSlices, kDim);
	for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
	int remainSliceFactor = sliceFactor;
	for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
		for (int i = 0; i < remainSliceFactor; ++i) {
			uint dimSliceFactor = remainSliceFactor - i;
			if (numWGs[kI] >= pow2(dimSliceFactor)) {
				numSlices[kI] = pow2(dimSliceFactor);
				remainSliceFactor -= dimSliceFactor;
				break;
			}
		}
	}

	UNIQUE_PTR_SIZE_T(sliceGlobalSize, kDim);
	for (int kI = 0; kI < kDim; ++kI) {
		sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
	}

	AtomicKernelSet aks = createAtomicKernelSet(requiredNumSlices);
	if (kDim == 1) {
		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
			EDKernel kernelSlice = createKernelCopy(k);
			kernelSlice->setGlobalSize(sliceGlobalSize.get());

			std::unique_ptr<size_t[]> sliceOffset = std::make_unique<size_t[]>(kDim);
			sliceOffset[0] = sliceGlobalSize[0] * dim0;
			kernelSlice->setOffset(sliceOffset.get());

			if (sliceFactor != 0 && dim0 == numSlices[0] - 1) {
				std::unique_ptr<size_t[]> remainGlobalSize = std::make_unique<size_t[]>(kDim);
				remainGlobalSize[0] = sliceGlobalSize[0];
				if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
					remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
				}
				kernelSlice->setGlobalSize(remainGlobalSize.get());
			}
			aks->addKernel(kernelSlice);
		}
	} else if (kDim == 2) {
		for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
			for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
				EDKernel kernelSlice = createKernelCopy(k);
				kernelSlice->setGlobalSize(sliceGlobalSize.get());
				std::unique_ptr<size_t[]> sliceOffset = std::make_unique<size_t[]>(kDim);
				sliceOffset[0] = sliceGlobalSize[0] * dim0;
				sliceOffset[1] = sliceGlobalSize[1] * dim1;
				kernelSlice->setOffset(sliceOffset.get());
				if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1)) { // if it's the edge slice
					std::unique_ptr<size_t[]> remainGlobalSize = std::make_unique<size_t[]>(kDim);
					remainGlobalSize[0] = sliceGlobalSize[0];
					remainGlobalSize[1] = sliceGlobalSize[1];
					if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
						remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
					}
					if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
						remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
					}
					kernelSlice->setGlobalSize(remainGlobalSize.get());
				}
				aks->addKernel(kernelSlice);
			}
		}
	} else {
		for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
			for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
				for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
					EDKernel kernelSlice = createKernelCopy(k);
					kernelSlice->setGlobalSize(sliceGlobalSize.get());
					std::unique_ptr<size_t[]> sliceOffset = std::make_unique<size_t[]>(kDim);
					sliceOffset[0] = sliceGlobalSize[0] * dim0;
					sliceOffset[1] = sliceGlobalSize[1] * dim1;
					sliceOffset[2] = sliceGlobalSize[2] * dim2;
					kernelSlice->setOffset(sliceOffset.get());

					if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
							|| dim2 == numSlices[2] - 1)) { // if it's the edge slice
						std::unique_ptr<size_t[]> remainGlobalSize = std::make_unique<size_t[]>(kDim);
						remainGlobalSize[0] = sliceGlobalSize[0];
						remainGlobalSize[1] = sliceGlobalSize[1];
						remainGlobalSize[2] = sliceGlobalSize[2];
						if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
							remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
						}
						if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
							remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
						}
						if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
							remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
						}
						kernelSlice->setGlobalSize(remainGlobalSize.get());
					}
					aks->addKernel(kernelSlice);
				}
			}
		}
	}
	return aks;
}

AtomicKernelSet EDCL::
slicingKernelPow2_extra(EDKernel k, uint sliceFactor) {
	if (k->localSize == nullptr) {
		PRINT_RED
		flushed_printf("Kernel %s local work size not set!\n", k->name.data());
		PRINT_DEFAULT
		return nullptr;
	}
	uint kDim = k->workDim; // kernel NDRange dimensions
	for (int kI = 0; kI < kDim; ++kI) {
		if (k->globalSize[kI] % k->localSize[kI] != 0) {
			PRINT_RED
			flushed_printf("Kernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n",
										 k->name.data(),
										 kI, k->globalSize[kI], k->localSize[kI]);
			PRINT_DEFAULT
			return nullptr;
		}
	}
	uint requiredNumSlices = pow2(sliceFactor);
	UNIQUE_PTR_SIZE_T(numWGs, kDim);
	uint maxNumSlices = 1;
	for (int kI = 0; kI < kDim; ++kI) {
		numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
		maxNumSlices *= numWGs[kI];
	}
	if (maxNumSlices < requiredNumSlices) {
		PRINT_RED
		flushed_printf("Slice factor %d exceeds the kernel global size. ", sliceFactor);
		flushed_printf("Number of WGs({");
		for (int kI = 0; kI < kDim; ++kI)
			flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
		flushed_printf("}):%d. Required number of slices: %d.\n", maxNumSlices, requiredNumSlices);
		PRINT_DEFAULT
		return nullptr;
	}

	UNIQUE_PTR_SIZE_T(numSlices, kDim);
	for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
	int remainSliceFactor = sliceFactor;
	for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
		for (int i = 0; i < remainSliceFactor; ++i) {
			uint dimSliceFactor = remainSliceFactor - i;
			if (numWGs[kI] >= pow2(dimSliceFactor)) {
				numSlices[kI] = pow2(dimSliceFactor);
				remainSliceFactor -= dimSliceFactor;
				break;
			}
		}
	}

	UNIQUE_PTR_SIZE_T(sliceGlobalSize, kDim);
	for (int kI = 0; kI < kDim; ++kI) {
		sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
	}

	AtomicKernelSet aks = createAtomicKernelSet(requiredNumSlices);
	if (kDim == 1) {
		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
			EDKernel kernelSlice = createKernelCopy(k);
			kernelSlice->setGlobalSize(sliceGlobalSize.get());
			UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
			sliceOffset[0] = sliceGlobalSize[0] * dim0;
			kernelSlice->setOffset(sliceOffset.get());
			aks->addKernel(kernelSlice);
			if (sliceFactor != 0) {
				if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
					EDKernel extraSlice = createKernelCopy(k);
					UNIQUE_PTR_SIZE_T(extraOffset, kDim);
					extraOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
					extraSlice->setOffset(extraOffset.get());
					UNIQUE_PTR_SIZE_T(extraGlobal, kDim);
					extraGlobal[0] = k->globalSize[0] - extraSlice->offset[0];
					extraSlice->setGlobalSize(extraGlobal.get());
					aks->addKernel(extraSlice);
				}
			}
		}
	} else if (kDim == 2) {
		for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
			for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
				EDKernel kernelSlice = createKernelCopy(k);
				kernelSlice->setGlobalSize(sliceGlobalSize.get());
				UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
				sliceOffset[0] = sliceGlobalSize[0] * dim0;
				sliceOffset[1] = sliceGlobalSize[1] * dim1;
				kernelSlice->setOffset(sliceOffset.get());
				aks->addKernel(kernelSlice);
				if (sliceFactor != 0) { // if it's the edge slice
					if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
						EDKernel extraSlice = createKernelCopy(k);
						UNIQUE_PTR_SIZE_T(remainGlobalSize, kDim);
						remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
						remainGlobalSize[1] = sliceGlobalSize[1];
						UNIQUE_PTR_SIZE_T(remainOffset, kDim);
						remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
						remainOffset[1] = kernelSlice->offset[1];
						extraSlice->setOffset(remainOffset.get());
						extraSlice->setGlobalSize(remainGlobalSize.get());
						aks->addKernel(extraSlice);
					}
					if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
						EDKernel extraSlice = createKernelCopy(k);
						UNIQUE_PTR_SIZE_T(remainGlobalSize, kDim);
						remainGlobalSize[0] = sliceGlobalSize[0];
						remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
						UNIQUE_PTR_SIZE_T(remainOffset, kDim);
						remainOffset[0] = kernelSlice->offset[0];
						remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
						extraSlice->setOffset(remainOffset.get());
						extraSlice->setGlobalSize(remainGlobalSize.get());
						aks->addKernel(extraSlice);
					}
					if ((dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1])
							&& dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
						EDKernel extraSlice = createKernelCopy(k);
						UNIQUE_PTR_SIZE_T(remainGlobalSize, kDim);
						remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
						remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
						UNIQUE_PTR_SIZE_T(remainOffset, kDim);
						remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
						remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];
						extraSlice->setOffset(remainOffset.get());
						extraSlice->setGlobalSize(remainGlobalSize.get());
						aks->addKernel(extraSlice);
					}
				}
			}
		}
	} else {
		for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
			for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
				for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
					EDKernel kernelSlice = createKernelCopy(k);
					kernelSlice->setGlobalSize(sliceGlobalSize.get());
					UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
					sliceOffset[0] = sliceGlobalSize[0] * dim0;
					sliceOffset[1] = sliceGlobalSize[1] * dim1;
					sliceGlobalSize[2] = sliceGlobalSize[2] * dim2;
					kernelSlice->setOffset(sliceOffset.get());
					if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
							|| dim2 == numSlices[2] - 1)) { // if it's the edge slice
						UNIQUE_PTR_SIZE_T(remainGlobalSize, kDim);
						remainGlobalSize[0] = sliceGlobalSize[0];
						remainGlobalSize[1] = sliceGlobalSize[1];
						remainGlobalSize[2] = sliceGlobalSize[2];
						if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
							remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
						}
						if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
							remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
						}
						if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
							remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
						}
						kernelSlice->setGlobalSize(remainGlobalSize.get());
					}
					aks->addKernel(kernelSlice);
				}
			}
		}
	}
	return aks;
}

AtomicKernelSet EDCL::
slicingKernelPow2_extra_equal(EDKernel k, uint sliceFactor) {
	if (k->localSize == nullptr) {
		flushed_printf("\033[1;31mKernel %s local work size not set!\n\033[0m", k->name.data());
		return nullptr;
	}
	uint kDim = k->workDim; // kernel NDRange dimensions
	for (int kI = 0; kI < kDim; ++kI) {
		if (k->globalSize[kI] % k->localSize[kI] != 0) {
			flushed_printf(
					"\033[1;31mKernel %s global work size is not divisible by local work size! Dim: %d, global:%d, local:%d\n\033[0m",
					k->name.data(),
					kI,
					k->globalSize[kI],
					k->localSize[kI]);
			return nullptr;
		}
	}

	uint requiredNumSlices = pow2(sliceFactor);
	auto *numWGs = new size_t[kDim];
	uint maxNumSlices = 1;
	for (int kI = 0; kI < kDim; ++kI) {
		numWGs[kI] = k->globalSize[kI] / k->localSize[kI];
		maxNumSlices *= numWGs[kI];
	}

//  for (int kI = 0; kI < k->workDim; ++kI) {
//	debug_printf(__func__, "%d: Kernel localSize[%d]: %d\n", __LINE__, kI, k->localSize[kI]);
//  }

	if (maxNumSlices < requiredNumSlices) {
		flushed_printf("\033[1;31mKernel %s slice factor %d exceeds the kernel global size. ",
									 k->name.c_str(),
									 sliceFactor);
		flushed_printf("Number of WGs({");
		for (int kI = 0; kI < kDim; ++kI)
			flushed_printf("%d%s", numWGs[kI], kI + 1 == kDim ? "" : ", ");
		flushed_printf("}):%d. Required number of slices: %d.\n\033[0m", maxNumSlices, requiredNumSlices);
		delete[] numWGs;
		return nullptr;
	}

	auto *numSlices = new size_t[kDim];
	for (int kI = 0; kI < kDim; ++kI) numSlices[kI] = 1;
	int remainSliceFactor = sliceFactor;
	for (int kI = 0; kI < kDim && remainSliceFactor > 0; ++kI) {
		for (int i = 0; i < remainSliceFactor; ++i) {
			uint dimSliceFactor = remainSliceFactor - i;
			if (numWGs[kI] >= pow2(dimSliceFactor)) {
				numSlices[kI] = pow2(dimSliceFactor);
				remainSliceFactor -= dimSliceFactor;
				break;
			}
		}
	}

	UNIQUE_PTR_SIZE_T(sliceGlobalSize, kDim);
	for (int kI = 0; kI < kDim; ++kI) {
		sliceGlobalSize[kI] = (numWGs[kI] / numSlices[kI]) * k->localSize[kI];
	}

	AtomicKernelSet aks = createAtomicKernelSet(requiredNumSlices);
	if (kDim == 1) {
		for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
			EDKernel kernelSlice = createKernelCopy(k);
			kernelSlice->setGlobalSize(sliceGlobalSize.get());
			UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
			sliceOffset[0] = sliceGlobalSize[0] * dim0;
			kernelSlice->setOffset(sliceOffset.get());
			aks->addKernel(kernelSlice);
			if (sliceFactor != 0) {
				size_t remainWI0 = k->globalSize[0] - (kernelSlice->offset[0] + sliceGlobalSize[0]);
				if (dim0 == numSlices[0] - 1 && remainWI0 > 0) {
					size_t remainChunk = ceil((float)remainWI0 / (float)sliceGlobalSize[0]);
					for (int i = 0; i < remainChunk; ++i) {
						EDKernel extraSlice = createKernelCopy(k);
						UNIQUE_PTR_SIZE_T(extraSliceOffset, kDim);
						extraSliceOffset[0] = kernelSlice->offset[0] + (i + 1) * sliceGlobalSize[0];
						extraSlice->setOffset(extraSliceOffset.get());
						UNIQUE_PTR_SIZE_T(extraSliceGlobal, kDim);
						extraSlice->setGlobalSize(extraSliceGlobal.get());
						extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
						aks->addKernel(extraSlice);
					}
				}
			}
		}
	} else if (kDim == 2) {
		for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
			for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
				EDKernel kernelSlice = createKernelCopy(k);
				kernelSlice->setGlobalSize(sliceGlobalSize.get());
				UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
				sliceOffset[0] = sliceGlobalSize[0] * dim0;
				sliceOffset[1] = sliceGlobalSize[1] * dim1;
//		debug_printf(__func__, "sliceOffset{%d, %d}\n", sliceOffset[0], sliceOffset[1]);
				kernelSlice->setOffset(sliceOffset.get());
//		debug_printf(__func__, "Kernel slice info ");
//		kernelSlice->printInfo();
				aks->addKernel(kernelSlice);
				if (sliceFactor != 0) { // if it's the edge slice
					if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
						size_t remainGlobalSize[2];
						remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
						remainGlobalSize[1] = sliceGlobalSize[1];
						size_t remainOffset[2];
						remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
						remainOffset[1] = kernelSlice->offset[1];

						size_t remainChunk = ceil((float)remainGlobalSize[0] / (float)sliceGlobalSize[0]);
						for (int i = 0; i < remainChunk; ++i) {
							EDKernel extraSlice = createKernelCopy(k);
							UNIQUE_PTR_SIZE_T(extraSliceOffset, kDim);
							extraSliceOffset[0] = remainOffset[0] + sliceGlobalSize[0] * i;
							extraSliceOffset[1] = remainOffset[1];
							extraSlice->setOffset(extraSliceOffset.get());
							UNIQUE_PTR_SIZE_T(extraSliceGlobal, kDim);
							extraSlice->setGlobalSize(extraSliceGlobal.get());
							extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
							extraSlice->globalSize[1] = remainGlobalSize[1];
							aks->addKernel(extraSlice);
						}

					}
					if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
						size_t remainGlobalSize[2];
						remainGlobalSize[0] = sliceGlobalSize[0];
						remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
						size_t remainOffset[2];
						remainOffset[0] = kernelSlice->offset[0];
						remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];

						size_t remainChunk = ceil((float)remainGlobalSize[1] / (float)sliceGlobalSize[1]);
						for (int i = 0; i < remainChunk; ++i) {
							EDKernel extraSlice = createKernelCopy(k);
							UNIQUE_PTR_SIZE_T(extraSliceOffset, kDim);
							extraSliceOffset[0] = remainOffset[0];
							extraSliceOffset[1] = remainOffset[1] + sliceGlobalSize[1] * i;
							extraSlice->setOffset(extraSliceOffset.get());
							UNIQUE_PTR_SIZE_T(extraSliceGlobal, kDim);
							extraSlice->setGlobalSize(extraSliceGlobal.get());
							extraSlice->globalSize[0] = remainGlobalSize[0];
							extraSlice->globalSize[1] = MIN(sliceGlobalSize[1], k->globalSize[1] - extraSlice->offset[1]);
							aks->addKernel(extraSlice);
						}
					}
					if ((dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1])
							&& dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
						size_t remainGlobalSize[2];
						remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0] - kernelSlice->globalSize[0];
						remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1] - kernelSlice->globalSize[1];
						size_t remainOffset[2];
						remainOffset[0] = kernelSlice->offset[0] + kernelSlice->globalSize[0];
						remainOffset[1] = kernelSlice->offset[1] + kernelSlice->globalSize[1];

						size_t remainChunk[2];
						remainChunk[0] = ceil((float)remainGlobalSize[0] / (float)sliceGlobalSize[0]);
						remainChunk[1] = ceil((float)remainGlobalSize[1] / (float)sliceGlobalSize[1]);
						for (int i = 0; i < remainChunk[0]; ++i) {
							for (int j = 0; j < remainChunk[1]; ++j) {
								EDKernel extraSlice = createKernelCopy(k);
								UNIQUE_PTR_SIZE_T(extraSliceOffset, kDim);
								extraSliceOffset[0] = remainOffset[0] + sliceGlobalSize[0] * i;
								extraSliceOffset[1] = remainOffset[1] + sliceGlobalSize[1] * j;
								UNIQUE_PTR_SIZE_T(extraSliceGlobal, kDim);
								extraSlice->setGlobalSize(extraSliceGlobal.get());
								extraSlice->globalSize[0] = MIN(sliceGlobalSize[0], k->globalSize[0] - extraSlice->offset[0]);
								extraSlice->globalSize[1] = MIN(sliceGlobalSize[1], k->globalSize[1] - extraSlice->offset[1]);
								aks->addKernel(extraSlice);
							}
						}
					}
				}
			}
		}
	} else {
		for (int dim2 = 0; dim2 < numSlices[2]; ++dim2) {
			for (int dim1 = 0; dim1 < numSlices[1]; ++dim1) {
				for (int dim0 = 0; dim0 < numSlices[0]; ++dim0) {
					EDKernel kernelSlice = createKernelCopy(k);
					kernelSlice->setGlobalSize(sliceGlobalSize.get());
					UNIQUE_PTR_SIZE_T(sliceOffset, kDim);
					sliceOffset[0] = sliceGlobalSize[0] * dim0;
					sliceOffset[1] = sliceGlobalSize[1] * dim1;
					sliceOffset[2] = sliceGlobalSize[2] * dim2;
					kernelSlice->setOffset(sliceOffset.get());
					if (sliceFactor != 0 && (dim0 == numSlices[0] - 1 || dim1 == numSlices[1] - 1
							|| dim2 == numSlices[2] - 1)) { // if it's the edge slice
						UNIQUE_PTR_SIZE_T(remainGlobalSize, kDim);
						remainGlobalSize[0] = sliceGlobalSize[0];
						remainGlobalSize[1] = sliceGlobalSize[1];
						remainGlobalSize[2] = sliceGlobalSize[2];
						if (dim0 == numSlices[0] - 1 && kernelSlice->offset[0] + sliceGlobalSize[0] < k->globalSize[0]) {
							remainGlobalSize[0] = k->globalSize[0] - kernelSlice->offset[0];
						}
						if (dim1 == numSlices[1] - 1 && kernelSlice->offset[1] + sliceGlobalSize[1] < k->globalSize[1]) {
							remainGlobalSize[1] = k->globalSize[1] - kernelSlice->offset[1];
						}
						if (dim2 == numSlices[2] - 1 && kernelSlice->offset[2] + sliceGlobalSize[2] < k->globalSize[2]) {
							remainGlobalSize[2] = k->globalSize[2] - kernelSlice->offset[2];
						}
						kernelSlice->setGlobalSize(remainGlobalSize.get());
					}
					aks->addKernel(kernelSlice);
				}
			}
		}
	}
	delete[] numWGs;
	delete[] numSlices;
	return aks;

}

__unused void EDCL::confirmAKSEnv(EDQueue Queue, const AtomicKernelSet &aks) const {
	for (auto &k:aks->kernels) {
		confirmExeEnv(Queue, k);
	}
}

