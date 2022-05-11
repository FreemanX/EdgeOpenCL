//
// Created by pfxu on 7/7/19.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_KERNELMANAGER_H
#define MOBILEHETEROGENOUSPROJECT_CLION_KERNELMANAGER_H

#include "heteroCompLib.h"

enum PARAM_TYPE {
    USERZBUFFER, FLOAT, INT, DOUBLE, LONG
};

typedef std::pair<PARAM_TYPE, void *> ParamPtr;
typedef std::map<int, ParamPtr> ParamMap;

struct HeteroKernel {
    cl_kernel *GPUKernel;

    void (*CPUKernel)(ParamMap *paramMap);
};

void *extractPtr(ParamMap *map, int pos);

class KernelManger {
public:
    KernelManger(cl_wrapper *cl);

    void executeKernelOnCPU();

    void executeKernelOnThread();

    void executeKernelOnGPU();

    void setHKernelArg(int pos, PARAM_TYPE type, void *paramPtr);

    HeteroKernel *getHeteroKernel();

    void setHeteroKernel(void (*CPUKernel)(ParamMap *paramMap), cl_kernel *GPUKernel);

    ~KernelManger();

private:
    cl_wrapper *cl;
    ParamMap paramMap;
    ParamMap cpuParamMap;
    HeteroKernel heteroKernel;
};

#endif //MOBILEHETEROGENOUSPROJECT_CLION_KERNELMANAGER_H
