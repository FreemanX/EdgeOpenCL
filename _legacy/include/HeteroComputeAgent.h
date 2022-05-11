//
// Created by pfxu on 7/5/19.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_HETEROCOMPUTEAGENT_H
#define MOBILEHETEROGENOUSPROJECT_CLION_HETEROCOMPUTEAGENT_H

#include "heteroCompLib.h"

class HeteroComputeAgent {

public:
    HeteroComputeAgent();

    ~HeteroComputeAgent();

    UserZBuffer applyZeroCopyBuffer(size_t size);


    void releaseZeroCopyBuffer(UserZBuffer *userZBuffer);

private:
    CLInfo clInfo{};
    cl_wrapper *cl{};
    MemoryManager *memoryManager;

    void initCLEnv();
};


#endif //MOBILEHETEROGENOUSPROJECT_CLION_HETEROCOMPUTEAGENT_H
