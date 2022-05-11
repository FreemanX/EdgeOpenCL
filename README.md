# Edge OpenCL Library  
EdgeOpenCL(EDCL) is a C++11 based framework developed to hide the complexity of programming on a heterogeneous mobile platform like Snapdragon 845 for developers.   
<img width="483" alt="EDCL" src="https://user-images.githubusercontent.com/9710644/167869887-29166cc9-47b6-47be-b399-b892f706f14c.png">   
The image above shows the structure of EDCL. EDCL holds two OpenCL runtime(RT) objects, one for the CPU and another for the GPU. Each OpenCL RT object keeps a dynamic library handler(returned by _dlopen()_) associated with corresponding OpenCL libray. An OpenCL function call layer is implemented to handle OpenCL function invocations. 
EDCL also keep track of the resource usage like creating and releasing command queues, buffer objects, kernel objects, etc.. Here's an example of using EDCL performing a vector addition followed by a matrix addition:
```C++
#include "EDCL.h"
#include "EDCL_Scheduler.h"

const char *source =
#include "vecAdd.cl"
unsigned int vDim = 4096;
uint n = vDim * vDim;
size_t nByte = n * sizeof(float);
size_t glob_MA[2] = {vDim / 4, vDim};
size_t local_MA[2] = {32, 32};
size_t glob_VA[1] = {n};
size_t local_VA[1] = {1024};

int main() {
  EDCL *edcl = new EDCL();
  // Init Buffers
  EDBuffer A = edcl->createBuffer(nByte);
  EDBuffer B = edcl->createBuffer(nByte);
  EDBuffer C = edcl->createBuffer(nByte);
  fillBuf<float>((float*)A->hostPtr, n);
  fillBuf<float>((float*)B->hostPtr, n);
  memset(C->hostPtr, 0, nByte);
  // Init Kernels
  EDProgram prog = edcl->createProgram(&source);
  // program, kernel name, number of args
  EDKernel matAdd = edcl->createKernel(prog, "matAdd", 5); 
  EDKernel vecAdd = edcl->createKernel(prog, "vecAdd", 4);
  // Config kernel: work_dim, global_work_size, local_work_size, arg1, arg2...
  matAdd->configKernel(2, glob_MA, local_MA, A, B, C, &vDim, &vDim);
  matAdd->addSubscriber(edcl, vecAdd);
  vecAdd->configKernel(1, glob_VA, local_VA, A, B, C, &n);
  // Get a scheduler instance
  auto scheduler = EDCL_Scheduler::getInstance(edcl);
  scheduler->submitKernels(matAdd, vecAdd);
  // Create scheduling strategy 
  Sequential strategy(edcl);
  // Set strategy parameters
  strategy.setExecutionDevice(GPU);
  // Execution
  scheduler->executeKernels(strategy);
  std::cout << "Planing overhead: " << strategy.planingTime;
  std::cout << ", execution time: " << strategy.executionTime << "\n";

  delete edcl;
  return 0;
}
```
Two kernels (vecAdd.cl is located in kernels/vecAdd.cl), _matAdd_ and _vecAdd_, are created and submitted to the scheduler. Here _vecAdd_ is a subscriber of _matAdd_ meaning that the execution of _vecAdd_ depends on the completion of _matAdd_. A simple sequential schedualing strategy is used here. Developers can implement their own scheduling strategies. 
## Main Components
### EDProgram and EDKernel
EDProgram and EDKernel both keep two copies of their corresponding OpenCL objects for GPU and CPU. EDKernel will keep a reference to each kernel argument and will only set the arguments once the execution environment is confirmed (with EDQueue). An EDKernel can accept other kernels as its subscribers and all the subscribers will wait for its completion of execution. For example, kernel A can subscribe to kernel B `B->addSubscriber(A)`, and A will wait for the completion of execution of Kernel B. To achieve cross-platform clevent like functionalities, an EDUserEvent which contains user events for both OpenCL platforms will be created. Later, kernel A will submit the corresponding user event to an event waiting list when calling _clEnqueueNDRangeKernel_. Once B completes, it will notify all its subscribers by setting the status of all EDUserEvents in the subscriber list to _COMPLETE_.    
### EDQueue   
EDQueue is a wrapper class for OpenCL's _cl\_command\_queue_. Each EDQueue is bound to a physical device. During the initialization stage, EDCL will gather basic information about the SoC like number of levels of heterogeneity, number of compute units(CUs) at each level. For a Snapdragon 845(the SoC I used), Qualcomm provides OpenCL library for the GPU. Though there are two CUs in Adreno 630, Qualcomm's OpenCL implementation doesn't allow users to create sub-devices. As for the CPU, you can use your own CPU OpenCL library on Android like _POCL_ that you can create sub-device containing any number of CUs and physically bind each CU to a CPU core by setting its affinity. As a result, various device partitioning scheme can be created and developers can assign different tasks to each device partition.    
![devicePartition copy](https://user-images.githubusercontent.com/9710644/167885279-d6e3f974-4e2d-4e76-bd3d-c182e959cf64.jpg)   
The image above demonstrates some possible partitioning examples which the framework supports. All user needs to do are creating EDQueue for OpenCL kernel execution and setting core affinities.
### EDBuffer
To take advantage of shared main memory design of the mobile SoC, we can choose to keep only one physical copy of data and let both CPU and GPU _cl\_mem_ object and a host pointer point to that memory region as shown in the figure below:   
<img width="491" alt="edbuffer" src="https://user-images.githubusercontent.com/9710644/167887915-df89d831-c3fc-4b91-ae25-ebac77315024.png">   
### Scheduler and Strategy
The scheduler class is responsible for managing kernels and executing those kernels with given strategy. The scheduling process consists of three phases: planing, executing and finishing. Developers only need to implement these three functions for their own strategies. An instance of resource manager will be passed to each phase so that a strategy can request system resources like a command queue for a specific device where no other kernels can be executed except via this command queue, setting frequencies of processors, etc.. During the planing phase, a strategy can prepare for its necessary data. For example, a machine learning based scheduling scheme may need to extract the kernel features and train some models. Runtime execution process should be implemented in the _execute()_ function. When the execution completes, resources such as temporary memory spaces and command queues should be released in the finishing phase.
### Resource Manager
This class provides some handy functions related to the hardware like tracking the execution, recording energy consumption indicator, setting CPU core affinity, setting processors' frequencies, querying hardware information, etc. 
## Kernel Execution
## Kernel Slicing
