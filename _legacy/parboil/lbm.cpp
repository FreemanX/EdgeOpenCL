//
// Created by pfxu on 20/8/2019.
//

#include "lbm.h"

int runCPU(const char *input, double *exeTime) {
    int t;
    MAIN_Param param{};
    MAIN_parseParam(input, &param);
    MAIN_printInfo(&param);
    CPU_MAIN_initialize(&param);
    double timer = getCurrentTime();
    for (t = 1; t <= param.nTimeSteps; t++) {
        if (param.simType == CHANNEL) {
            LBM_handleInOutFlow(*srcGrid);
        }
        LBM_performStreamCollide(*srcGrid, *dstGrid);
        LBM_swapGrids(&srcGrid, &dstGrid);
    }
    timer = getCurrentTime() - timer;
    flushed_printf("\tlbm CPU done, time: %f sec\n", timer);
    *exeTime = timer;

    LBM_freeGrid((float **) &srcGrid);
    LBM_freeGrid((float **) &dstGrid);
    return 0;
}

int runGPU(const char *input, double *exeTime) {
    int t;
    MAIN_Param param{};
    MAIN_parseParam(input, &param);
    ZeroCopyMem<float> srcGrid_z{}, dstGrid_z{};
    MAIN_printInfo(&param);
    initCL();
    GLBM_allocateGrid(&srcGrid_z);
    GLBM_allocateGrid(&dstGrid_z);
    LBM_initializeGrid(srcGrid_z.hostPtr);
    LBM_initializeGrid(dstGrid_z.hostPtr);
    if (param.obstacleFilename != nullptr) {
        LBM_loadObstacleFile(srcGrid_z.hostPtr, param.obstacleFilename);
        LBM_loadObstacleFile(dstGrid_z.hostPtr, param.obstacleFilename);
    }
    LBM_initializeSpecialCellsForLDC(srcGrid_z.hostPtr);
    LBM_initializeSpecialCellsForLDC(dstGrid_z.hostPtr);
    LBM_showGridStatistics(srcGrid_z.hostPtr);
    double timer = getCurrentTime();

    for (t = 1; t <= param.nTimeSteps; t++) {
        OpenCL_LBM_performStreamCollide(&srcGrid_z.deviceBuffer, &dstGrid_z.deviceBuffer);
        LBM_swapGrids(&srcGrid_z, &dstGrid_z);
    }

    timer = getCurrentTime() - timer;
    flushed_printf("\tlbm GPU done, time: %f sec\n", timer);
    *exeTime = timer;

    release_zero_copy_region(&srcGrid_z);
    release_zero_copy_region(&dstGrid_z);
    return 0;
}

int main(int argc, char **argv) {
    char const *input = "datasets/lbm/long/input/120_120_150_ldc.of";
    char const *input_small = "datasets/lbm/short/input/120_120_150_ldc.of";
    benchmark(input, input_small, runCPU, runGPU);

    return 0;
}
