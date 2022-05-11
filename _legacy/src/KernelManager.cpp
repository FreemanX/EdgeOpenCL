//
// Created by pfxu on 7/7/19.
//

#include "KernelManager.h"
#include <map>


KernelManger::KernelManger(cl_wrapper *cl) {
    this->cl = cl;
}

KernelManger::~KernelManger() {
    paramMap.clear();
}

void KernelManger::setHKernelArg(int pos, PARAM_TYPE type, void *paramPtr) {
    ParamPtr hPtr = ParamPtr(type, paramPtr);
    paramMap[pos] = hPtr;
//    auto const hResult = paramMap.insert(std::pair<int, ParamPtr>(pos, hPtr));
//    if (!hResult.second) hResult.first->second = hPtr;
    if (type == USERZBUFFER) {

    }
}

void KernelManger::setHeteroKernel(void (*CPUKernel)(ParamMap *), cl_kernel *GPUKernel) {
    heteroKernel.CPUKernel = CPUKernel;
    heteroKernel.GPUKernel = GPUKernel;
}

HeteroKernel *KernelManger::getHeteroKernel() {
    return &heteroKernel;
}

void KernelManger::executeKernelOnCPU() {
    heteroKernel.CPUKernel(&cpuParamMap);
}

