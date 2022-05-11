//
// Created by pfxu on 7/5/19.
//
#include "HeteroComputeAgent.h"

void HeteroComputeAgent::initCLEnv() {
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices->deviceId);
}

HeteroComputeAgent::HeteroComputeAgent() {
    initCLEnv();
    memoryManager = new MemoryManager(cl);
}

HeteroComputeAgent::~HeteroComputeAgent() {
}

UserZBuffer HeteroComputeAgent::applyZeroCopyBuffer(size_t size) {
    return memoryManager->applyZeroCopyBuffer(size);
}

void HeteroComputeAgent::releaseZeroCopyBuffer(UserZBuffer *userZBuffer) {
    memoryManager->releaseUserZBuffer(userZBuffer);
}
