//
// Created by pfxu on 7/19/19.
//
#include "benchmark.h"
#include "mri-gridding.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//=========================================main=====================================//

int runCPU(const char *input, double *exeTime) {
    char uksdata[250];
    parameters params;

    FILE *uksfile_f = nullptr;
    FILE *uksdata_f = nullptr;

    strcpy(uksdata, input);
    strcat(uksdata, ".data");

    uksfile_f = fopen(input, "r");
    if (uksfile_f == nullptr) {
        printf("ERROR: Could not open %s\n", input);
        return 1;
    }
    printf("\nReading parameters\n");
    params.binsize = 256;
    setParameters(uksfile_f, &params);
    auto *samples = (ReconstructionSample *) malloc(params.numSamples * sizeof(ReconstructionSample)); //Input Data
    float *LUT; //use look-up table for faster execution on CPU (intermediate data)
    unsigned int sizeLUT; //set in the function calculateLUT (intermediate data)

    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];

    auto *gridData = (cmplx *) calloc(gridNumElems, sizeof(cmplx)); //Output Data
    auto *sampleDensity = (float *) calloc(gridNumElems, sizeof(float)); //Output Data

    if (samples == nullptr) {
        printf("ERROR: Unable to allocate memory for input data\n");
        return 1;
    }
    if (sampleDensity == nullptr || gridData == nullptr) {
        printf("ERROR: Unable to allocate memory for output data\n");
        return 1;
    }
    uksdata_f = fopen(uksdata, "rb");
    if (uksdata_f == nullptr) {
        printf("ERROR: Could not open data file\n");
        return 1;
    }
    printf("Reading input data from files\n");
    unsigned int n = readSampleData(params, uksdata_f, samples);
    fclose(uksdata_f);
    if (params.useLUT) {
        printf("Generating Look-Up Table\n");
        float beta = PI * sqrt(4 * params.kernelWidth * params.kernelWidth / (params.oversample * params.oversample) *
                               (params.oversample - .5) * (params.oversample - .5) - .8);
        calculateLUT(beta, params.kernelWidth, &LUT, &sizeLUT);
    }

    flushed_printf("\n\tmri-gridding Running...\n");
    double timer = getCurrentTime();
    gridding_Gold(n, params, samples, LUT, sizeLUT, gridData, sampleDensity);
    timer = getCurrentTime() - timer;
    flushed_printf("\tmri-gridding CPU done, time: %f sec\n", timer);

    *exeTime = timer;
    if (params.useLUT) {
        free(LUT);
    }
    free(samples);
    free(gridData);
    free(sampleDensity);

    return 0;
}


static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_kernel binning_kernel;
cl_kernel reorder_kernel;
cl_kernel gridding_GPU;
cl_kernel splitSort;
cl_kernel splitRearrange;
cl_kernel scan_L1_kernel;
cl_kernel scan_inter1_kernel;
cl_kernel scan_inter2_kernel;
cl_kernel uniformAdd;
cl_command_queue exeCmdQueue;

void initCL() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    exeCmdQueue = cl->createProfilingCmdQueue();
}

#define GRID_SIZE 65535
#define NUM_BANKS 16
#define LOG_NUM_BANKS 4
#define EXPANDED_SIZE(__x) (__x+(__x>>LOG_NUM_BANKS)+(__x>>(2*LOG_NUM_BANKS)))

void scanLargeArray(unsigned int gridNumElems, cl_mem *data_d, size_t *workItemSizes) {
    size_t blockSize = (workItemSizes[0] * 2 < 1024) ? workItemSizes[0] * 2 : 1024;
    // Run the prescan
    unsigned int size = (gridNumElems + blockSize - 1) / blockSize;
    unsigned int dim_block = 0;
    unsigned int current_max = size * blockSize;
    for (int block_size_lcv = 128; block_size_lcv <= blockSize; block_size_lcv *= 2) {
        unsigned int array_size = block_size_lcv;
        while (array_size < size) {
            array_size *= block_size_lcv;
        }
        if (array_size <= current_max) {
            current_max = array_size;
            dim_block = block_size_lcv;
        }
    }
    cl_int ciErrNum;
    cl_mem inter_d;
    unsigned int *zeroData;
    zeroData = (unsigned int *) calloc(current_max, sizeof(unsigned int));
    if (zeroData == nullptr) {
        fprintf(stderr, "Could not allocate host memory! (%s)\n", __FILE__);
        exit(1);
    }
    inter_d = clCreateBuffer(cl->context, CL_MEM_COPY_HOST_PTR,
                             current_max * sizeof(unsigned int), zeroData, &ciErrNum);
    cl_wrapper::checkError(ciErrNum);

    free(zeroData);

    OCL_ERRCK_RETVAL(clSetKernelArg(scan_L1_kernel, 1, sizeof(cl_mem), (void *) data_d));
    OCL_ERRCK_RETVAL(clSetKernelArg(scan_L1_kernel, 3, sizeof(cl_mem), (void *) &inter_d));

    OCL_ERRCK_RETVAL(clSetKernelArg(scan_inter1_kernel, 0, sizeof(cl_mem), (void *) &inter_d));
    OCL_ERRCK_RETVAL(clSetKernelArg(scan_inter2_kernel, 0, sizeof(cl_mem), (void *) &inter_d));

    OCL_ERRCK_RETVAL(clSetKernelArg(uniformAdd, 1, sizeof(cl_mem), (void *) data_d));
    OCL_ERRCK_RETVAL(clSetKernelArg(uniformAdd, 3, sizeof(cl_mem), (void *) &inter_d));

    cl_event event;
    for (unsigned int i = 0; i < (size + GRID_SIZE - 1) / GRID_SIZE; i++) {
        unsigned int gridSize = ((size - (i * GRID_SIZE)) > GRID_SIZE) ? GRID_SIZE : (size - i * GRID_SIZE);
        unsigned int numElems = ((gridNumElems - (i * GRID_SIZE * blockSize)) > (GRID_SIZE * blockSize)) ?
                                (GRID_SIZE * blockSize) : (gridNumElems - (i * GRID_SIZE * blockSize));
        unsigned int data_offset = i * GRID_SIZE * blockSize;
        unsigned int inter_offset = i * GRID_SIZE;
        OCL_ERRCK_RETVAL(clSetKernelArg(scan_L1_kernel, 0, sizeof(unsigned int), &numElems));
        OCL_ERRCK_RETVAL(clSetKernelArg(scan_L1_kernel, 2, sizeof(unsigned int), &data_offset));
        OCL_ERRCK_RETVAL(clSetKernelArg(scan_L1_kernel, 4, sizeof(unsigned int), &inter_offset));
        size_t block[1] = {blockSize / 2};
        size_t grid[1] = {gridSize * block[0]};
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, scan_L1_kernel, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);
    }

    unsigned int stride = 1;
    for (unsigned int d = current_max; d > 1; d /= dim_block) {
        size_t block[1] = {dim_block / 2};
        size_t grid[1] = {(d / dim_block) * block[0]};
        OCL_ERRCK_RETVAL(clSetKernelArg(scan_inter1_kernel, 1, sizeof(unsigned int), &stride));
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, scan_inter1_kernel, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);
        stride *= dim_block;
    }
    unsigned int singleZero = 0;
    OCL_ERRCK_RETVAL(clEnqueueWriteBuffer(exeCmdQueue, inter_d, CL_TRUE,
                                          (current_max - 1) * sizeof(unsigned int), // Offset in bytes
                                          sizeof(unsigned int), // Size of data to write
                                          &singleZero, // Host Source
                                          0, nullptr, &event));
    clWaitForEvents(1, &event);
    for (unsigned int d = dim_block; d <= current_max; d *= dim_block) {
        stride /= dim_block;
        size_t block[1] = {dim_block / 2};
        size_t grid[1] = {(d / dim_block) * block[0]};
        OCL_ERRCK_RETVAL(clSetKernelArg(scan_inter2_kernel, 1, sizeof(unsigned int), &stride));
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, scan_inter2_kernel, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);
    }

    for (unsigned int i = 0; i < (size + GRID_SIZE - 1) / GRID_SIZE; i++) {
        unsigned int gridSize = ((size - (i * GRID_SIZE)) > GRID_SIZE) ? GRID_SIZE : (size - i * GRID_SIZE);
        unsigned int numElems = ((gridNumElems - (i * GRID_SIZE * blockSize)) > (GRID_SIZE * blockSize)) ?
                                (GRID_SIZE * blockSize) : (gridNumElems - (i * GRID_SIZE * blockSize));
        unsigned int data_offset = i * GRID_SIZE * blockSize;
        unsigned int inter_offset = i * GRID_SIZE;
        OCL_ERRCK_RETVAL(clSetKernelArg(uniformAdd, 0, sizeof(unsigned int), &numElems));
        OCL_ERRCK_RETVAL(clSetKernelArg(uniformAdd, 2, sizeof(unsigned int), &data_offset));
        OCL_ERRCK_RETVAL(clSetKernelArg(uniformAdd, 4, sizeof(unsigned int), &inter_offset));
        size_t block[1] = {blockSize / 2};
        size_t grid[1] = {gridSize * block[0]};
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, uniformAdd, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);
    }
    clReleaseEvent(event);
    OCL_ERRCK_RETVAL(clReleaseMemObject(inter_d));
}

#define UINT32_MAX 4294967295
#define BITS 4
#define LNB 4
#define SORT_BS 256

void sort(int numElems, unsigned int max_value, cl_mem *dkeysPtr, cl_mem *dvaluesPtr, cl_mem *dkeys_oPtr,
          cl_mem *dvalues_oPtr, size_t *workItemSizes) {
    size_t block[1] = {SORT_BS};
    size_t grid[1] = {((numElems + 4 * SORT_BS - 1) / (4 * SORT_BS)) * block[0]};

    unsigned int iterations = 0;
    while (max_value > 0) {
        max_value >>= BITS;
        iterations++;
    }

    unsigned int *zeroData;
    zeroData = (unsigned int *) calloc((1 << BITS) * grid[0], sizeof(unsigned int));
    if (zeroData == nullptr) {
        fprintf(stderr, "Could not allocate host memory! (%s: %d)\n", __FILE__, __LINE__);
        exit(1);
    }

    ZeroCopyMem<unsigned int> histo_z = init_zero_copy_region<unsigned int>
            (cl, (1 << BITS) * ((numElems + 4 * SORT_BS - 1) / (4 * SORT_BS)));
    memcpy(histo_z.hostPtr, zeroData,
           (1 << BITS) * ((numElems + 4 * SORT_BS - 1) / (4 * SORT_BS)) * sizeof(unsigned int));
    free(zeroData);

    OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 0, sizeof(int), &numElems));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 2, sizeof(cl_mem), (void *) dkeysPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 3, sizeof(cl_mem), (void *) dvaluesPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 4, sizeof(cl_mem), (void *) &histo_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 0, sizeof(int), &numElems));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 2, sizeof(cl_mem), (void *) dkeysPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 3, sizeof(cl_mem), (void *) dkeys_oPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 4, sizeof(cl_mem), (void *) dvaluesPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 5, sizeof(cl_mem), (void *) dvalues_oPtr));
    OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 6, sizeof(cl_mem), (void *) &histo_z.deviceBuffer));

    cl_event event;
    for (int i = 0; i < iterations; i++) {
        OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 1, sizeof(int), &i));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 2, sizeof(cl_mem), (void *) dkeysPtr));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitSort, 3, sizeof(cl_mem), (void *) dvaluesPtr));
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, splitSort, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);

        scanLargeArray(((numElems + 4 * SORT_BS - 1) / (4 * SORT_BS)) * (1 << BITS),
                       &histo_z.deviceBuffer, workItemSizes);

        OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 1, sizeof(int), &i));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 2, sizeof(cl_mem), (void *) dkeysPtr));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 3, sizeof(cl_mem), (void *) dkeys_oPtr));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 4, sizeof(cl_mem), (void *) dvaluesPtr));
        OCL_ERRCK_RETVAL(clSetKernelArg(splitRearrange, 5, sizeof(cl_mem), (void *) dvalues_oPtr));
        OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, splitRearrange, 1, 0, grid, block, 0, 0, &event));
        clWaitForEvents(1, &event);

//        cl_mem *temp = dkeysPtr;
//        dkeysPtr = dkeys_oPtr;
//        dkeys_oPtr = temp;
//
//        temp = dvaluesPtr;
//        dvaluesPtr = dvalues_oPtr;
//        dvalues_oPtr = temp;
    }
    OCL_ERRCK_RETVAL(clReleaseMemObject(histo_z.deviceBuffer));
}

int runGPU(const char *input, double *exeTime) {
    initCL();
    char uksdata[250];
    parameters params;

    FILE *uksfile_f = nullptr;
    FILE *uksdata_f = nullptr;

    strcpy(uksdata, input);
    strcat(uksdata, ".data");

    uksfile_f = fopen(input, "r");
    if (uksfile_f == nullptr) {
        printf("ERROR: Could not open %s\n", input);
        return 1;
    }
    printf("\nReading parameters\n");
    params.binsize = 256;
    setParameters(uksfile_f, &params);

    float *LUT = nullptr; //use look-up table for faster execution on CPU (intermediate data)
    unsigned int sizeLUT; //set in the function calculateLUT (intermediate data)
//    cmplx *gridData_gold; //Gold Output Data
//    float *sampleDensity_gold; //Gold Output Data

    size_t max_alloc_size = 0;
    (void) clGetDeviceInfo(cl->device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(size_t), &max_alloc_size, 0);
    size_t global_mem_size = 0;
    (void) clGetDeviceInfo(cl->device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(size_t), &global_mem_size, 0);
    size_t samples_size = params.numSamples * sizeof(ReconstructionSample);
    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];
    size_t output_size = gridNumElems * sizeof(cmplx);

    if (((samples_size + output_size) > global_mem_size) ||
        (samples_size > max_alloc_size) ||
        (output_size > max_alloc_size)) {
        fprintf(stderr, "Memory requirements for this dataset exceed device capabilities\n");
        exit(1);
    }

    cl_uint workItemDimensions;
    OCL_ERRCK_RETVAL(clGetDeviceInfo(
            cl->device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &workItemDimensions, nullptr));
    size_t workItemSizes[workItemDimensions];
    OCL_ERRCK_RETVAL(clGetDeviceInfo(
            cl->device, CL_DEVICE_MAX_WORK_ITEM_SIZES, workItemDimensions * sizeof(size_t), workItemSizes, nullptr));

    uksdata_f = fopen(uksdata, "rb");
    if (uksdata_f == NULL) {
        printf("ERROR: Could not open data file\n");
        exit(1);
    }

    printf("\tReading input data from files\n");

    ZeroCopyMem<ReconstructionSample> samples_z = init_zero_copy_region<ReconstructionSample>(cl, params.numSamples);
    unsigned int n = readSampleData(params, uksdata_f, samples_z.hostPtr);
    fclose(uksdata_f);

    if (params.useLUT) {
        printf("Generating Look-Up Table\n");
        float beta = PI * sqrt(4 * params.kernelWidth * params.kernelWidth / (params.oversample * params.oversample) *
                               (params.oversample - .5) * (params.oversample - .5) - .8);
        calculateLUT(beta, params.kernelWidth, &LUT, &sizeLUT);
    }

//    gridData_gold = (cmplx *) calloc(gridNumElems, sizeof(cmplx));
//    sampleDensity_gold = (float *) calloc(gridNumElems, sizeof(float));
//    if (sampleDensity_gold == NULL || gridData_gold == NULL) {
//        printf("ERROR: Unable to allocate memory for output data\n");
//        exit(1);
//    }

//    printf("\trunning gold version\n");
//
//    omp_set_num_threads(4);
//    std::vector<int> bigcpulist;
//    bigcpulist.reserve(4);
//    for (int i = 0; i < 4; ++i) {
//        bigcpulist.push_back(i + 4);
//    }
//    setcurthreadaffinity(bigcpulist);
//    gridding_Gold(n, params, samples_z.hostPtr, LUT, sizeLUT, gridData_gold, sampleDensity_gold);

    printf("\tRunning OpenCL version\n");

    ZeroCopyMem<cmplx> gridData_z = init_zero_copy_region<cmplx>(cl, gridNumElems);
    ZeroCopyMem<float> sampleDensity_z = init_zero_copy_region<float>(cl, gridNumElems);

    size_t blockSize = workItemSizes[0];
    int dims[3] = {8, 4, 2}; //size of a gridding block on the GPU

    /* x, y, z dimensions of the output grid (gridData) */
    int size_x = params.gridSize[0];
    int size_y = params.gridSize[1];
    int size_z = params.gridSize[2];
    int size_xy = size_y * size_x;

    gridNumElems = size_x * size_y * size_z;  // Total number of grid points
    float beta = PI * sqrt(4 * params.kernelWidth * params.kernelWidth / (params.oversample * params.oversample) *
                           (params.oversample - .5) * (params.oversample - .5) - .8);
    float cutoff = float(params.kernelWidth) / 2.0; // cutoff radius
    float cutoff2 = cutoff * cutoff;                // square of cutoff radius
    float _1overCutoff2 = 1 / cutoff2;              // 1 over square of cutoff radius

    unsigned int *zeroData = NULL, *maxIntData = NULL;
    size_t sizeZeroData = sizeof(float) * 2 * gridNumElems;
    if (n * sizeof(ReconstructionSample) > sizeZeroData)
        sizeZeroData = n * sizeof(ReconstructionSample);
    if ((sizeof(unsigned int) * (gridNumElems + 1)) > sizeZeroData)
        sizeZeroData = sizeof(unsigned int) * (gridNumElems + 1);
    if ((((n + 3) / 4) * 4) * sizeof(unsigned int) > sizeZeroData)
        sizeZeroData = (((n + 3) / 4) * 4) * sizeof(unsigned int);

    zeroData = (unsigned int *) malloc(sizeZeroData);
    maxIntData = (unsigned int *) malloc((((n + 3) / 4) * 4) * sizeof(unsigned int));
    memset(zeroData, 0, sizeZeroData);
    memset(maxIntData + n, 0xFF, (((n + 3) & ~(3)) - n) * sizeof(unsigned int));

    ZeroCopyMem<ReconstructionSample> sortedSample_z = init_zero_copy_region<ReconstructionSample>(cl, n);
    memcpy(sortedSample_z.hostPtr, zeroData, n * sizeof(ReconstructionSample));

    ZeroCopyMem<unsigned int> binStartAddr_z =
            init_zero_copy_region<unsigned int>(cl, gridNumElems + 1);
    memcpy(binStartAddr_z.hostPtr, zeroData, gridNumElems + 1 * sizeof(unsigned int));

    ZeroCopyMem<unsigned int> idxKey_z = init_zero_copy_region<unsigned int>(cl, ((n + 3) / 4) * 4);
    memcpy(idxKey_z.hostPtr, maxIntData, ((n + 3) / 4) * 4 * sizeof(unsigned int));

    ZeroCopyMem<unsigned int> idxValue_z = init_zero_copy_region<unsigned int>(cl, ((n + 3) / 4) * 4);
    memcpy(idxValue_z.hostPtr, maxIntData, ((n + 3) / 4) * 4 * sizeof(unsigned int));

    size_t blockSize_tmp = (workItemSizes[0] * 2 < 1024) ? workItemSizes[0] * 2 : 1024;
    // Run the prescan
    unsigned int size = (gridNumElems + blockSize - 1) / blockSize;

    unsigned int dim_block = 0;
    unsigned int current_max = size * blockSize_tmp;
    for (int block_size_lcv = 128; block_size_lcv <= blockSize_tmp; block_size_lcv *= 2) {
        unsigned int array_size = block_size_lcv;
        while (array_size < size) {
            array_size *= block_size_lcv;
        }
        if (array_size <= current_max) {
            current_max = array_size;
            dim_block = block_size_lcv;
        }
    }
    char compileOptions[1024];
    sprintf(compileOptions,
            "-D DYN_LOCAL_MEM_SIZE=%lu -D CUTOFF2_VAL=%f -D CUTOFF_VAL=%f -D CEIL_CUTOFF_VAL=%f -D GRIDSIZE_VAL1=%d -D GRIDSIZE_VAL2=%d -D GRIDSIZE_VAL3=%d -D SIZE_XY_VAL=%d -D ONE_OVER_CUTOFF2_VAL=%f",
            EXPANDED_SIZE(dim_block) * sizeof(unsigned int),
            cutoff2, cutoff, ceil(cutoff),
            params.gridSize[0], params.gridSize[1], params.gridSize[2],
            size_xy, _1overCutoff2);
    program = cl->createProgramWithOptions(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN, compileOptions);
    binning_kernel = cl->createKernel(program, "binning_kernel");
    reorder_kernel = cl->createKernel(program, "reorder_kernel");
    gridding_GPU = cl->createKernel(program, "gridding_GPU");
    splitSort = cl->createKernel(program, "splitSort");
    splitRearrange = cl->createKernel(program, "splitRearrange");
    scan_L1_kernel = cl->createKernel(program, "scan_L1_kernel");
    scan_inter1_kernel = cl->createKernel(program, "scan_inter1_kernel");
    scan_inter2_kernel = cl->createKernel(program, "scan_inter2_kernel");
    uniformAdd = cl->createKernel(program, "uniformAdd");
    free(maxIntData);

    size_t block1[1] = {blockSize};
    size_t grid1[1] = {((n + blockSize - 1) / blockSize) * block1[0]};
    cl_event event;

    double timer = getCurrentTime();
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 0, sizeof(unsigned int), &n));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 1, sizeof(cl_mem), (void *) &samples_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 2, sizeof(cl_mem), (void *) &idxKey_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 3, sizeof(cl_mem), (void *) &idxValue_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 4, sizeof(cl_mem), (void *) &binStartAddr_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 5, sizeof(int), &(params.binsize)));
    OCL_ERRCK_RETVAL(clSetKernelArg(binning_kernel, 6, sizeof(unsigned int), &gridNumElems));
    OCL_ERRCK_RETVAL(clSetKernelArg(reorder_kernel, 0, sizeof(unsigned int), &n));
    OCL_ERRCK_RETVAL(clSetKernelArg(reorder_kernel, 2, sizeof(cl_mem), (void *) &samples_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(reorder_kernel, 3, sizeof(cl_mem), (void *) &sortedSample_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, binning_kernel, 1, 0, grid1, block1, 0, 0, &event));
    clWaitForEvents(1, &event);
    ZeroCopyMem<unsigned int> keys_oz = init_zero_copy_region<unsigned int>(cl, n);
    ZeroCopyMem<unsigned int> values_oz = init_zero_copy_region<unsigned int>(cl, n);
    sort(n, gridNumElems + 1, &idxKey_z.deviceBuffer, &idxValue_z.deviceBuffer,
         &keys_oz.deviceBuffer, &values_oz.deviceBuffer, workItemSizes);
    clReleaseMemObject(idxKey_z.deviceBuffer);
    clReleaseMemObject(idxValue_z.deviceBuffer);
    OCL_ERRCK_RETVAL(clSetKernelArg(reorder_kernel, 1, sizeof(cl_mem), (void *) &values_oz.deviceBuffer));
    OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, reorder_kernel, 1, 0, grid1, block1, 0, 0, &event));
    clWaitForEvents(1, &event);
    OCL_ERRCK_RETVAL(clReleaseMemObject(keys_oz.deviceBuffer));
    OCL_ERRCK_RETVAL(clReleaseMemObject(values_oz.deviceBuffer));
    scanLargeArray(gridNumElems + 1, &binStartAddr_z.deviceBuffer, workItemSizes);

    unsigned int cpuStart = binStartAddr_z.hostPtr[gridNumElems];
    int CPUbin_size = int(n) - int(cpuStart);

    ReconstructionSample *CPUbin = &sortedSample_z.hostPtr[cpuStart];
    free(zeroData);

    size_t block2[3] = {static_cast<size_t>(dims[0]), static_cast<size_t>(dims[1]), static_cast<size_t>(dims[2])};
    size_t grid2[3] = {(size_x / dims[0]) * block2[0],
                       ((size_y * size_z) / (dims[1] * dims[2])) * block2[1], 1 * block2[2]};

    OCL_ERRCK_RETVAL(clSetKernelArg(gridding_GPU, 0, sizeof(cl_mem), (void *) &sortedSample_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(gridding_GPU, 1, sizeof(cl_mem), (void *) &binStartAddr_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(gridding_GPU, 2, sizeof(cl_mem), (void *) &gridData_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(gridding_GPU, 3, sizeof(cl_mem), (void *) &sampleDensity_z.deviceBuffer));
    OCL_ERRCK_RETVAL(clSetKernelArg(gridding_GPU, 4, sizeof(float), &beta));

    OCL_ERRCK_RETVAL(clEnqueueNDRangeKernel(exeCmdQueue, gridding_GPU, 3, 0, grid2, block2, 0, 0, &event));
    clWaitForEvents(1, &event);
    OCL_ERRCK_RETVAL(clReleaseMemObject(binStartAddr_z.deviceBuffer));

    gridding_Gold(CPUbin_size, params, CPUbin, LUT, sizeLUT, gridData_z.hostPtr, sampleDensity_z.hostPtr);
//    printf("\tDensity diff: %f\n", calVecDiff(sampleDensity_gold, sampleDensity_z.hostPtr, gridNumElems));

    timer = getCurrentTime() - timer;
    flushed_printf("\tmri-gridding GPU done, time: %f sec\n", timer);
    *exeTime = timer;

    if (params.useLUT)
        free(LUT);
    clReleaseMemObject(samples_z.deviceBuffer);
    clReleaseMemObject(gridData_z.deviceBuffer);
    clReleaseMemObject(sampleDensity_z.deviceBuffer);
//    free(gridData_gold);
//    free(sampleDensity_gold);


    return 0;
}

int main(int argc, char **argv) {
    char const *input = "datasets/mri-gridding/small/input/small.uks";
//    char const *input_small = "datasets/mri-gridding/mrig_small/small/input/small.uks";
//    double test;
//    runGPU(input, &test);
    benchmark(input, nullptr, runCPU, runGPU);

    return 0;
}
