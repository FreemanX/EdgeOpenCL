//
// Created by pfxu on 28/8/2019.
//
#include "HeteroComputeAgent.h"
#include <sys/mman.h>

double zeroCopyBench(cl_wrapper *cl, size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    ZeroCopyMem<float> tmp_a = init_zero_copy_region_byte<float>(cl, dataLength);
    time = getCurrentTime() - time;
    memset(&tmp_a.hostPtr[0], 1, dataLength);
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a.hostPtr[i] = 1;
    clReleaseMemObject(tmp_a.deviceBuffer);
    return time;
//    return 0;
}

double mallocBench(size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    auto *tmp_a = (float *) malloc(dataLength);
    time = getCurrentTime() - time;
    memset(&tmp_a[0], 1, dataLength);
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a[i] = 1;
    free(tmp_a);
    return time;
//    return 0;
}

double callocBench(size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    auto *tmp_a = (float *) calloc(dataLength, sizeof(float));
    time = getCurrentTime() - time;
    memset(&tmp_a[0], 0, dataLength * sizeof(float));
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a[i] = 1;
    free(tmp_a);
    return time;
    return 0;
}

double mmapBench(size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    char *ptr = reinterpret_cast<char *>(0x8000000);
//    mmap_allocate(ptr, dataLength / PAGE_SIZE);

    char *ret = (char *) mmap((char *) ptr, dataLength, PROT_EXEC | PROT_READ | PROT_WRITE,
                              MAP_NORESERVE | MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS, 0, 0);
    if (ret == MAP_FAILED) {
        fprintf(stderr, "Error mmap. errno: %d\n", errno);
        exit(-1);
    }

    time = getCurrentTime() - time;

    memset(&ptr[0], 1, dataLength);
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a[i] = 1;
    munmap(ptr, dataLength);
    return time;
    return 0;
}


double zeroCopyVecBench(cl_wrapper *cl, size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    ZeroCopyMem<float> tmp_a = init_zero_copy_region_byte<float>(cl, dataLength);
    ZeroCopyMem<float> tmp_b = init_zero_copy_region_byte<float>(cl, dataLength);
    ZeroCopyMem<float> tmp_c = init_zero_copy_region_byte<float>(cl, dataLength);
    time = getCurrentTime() - time;
    memset(&tmp_a.hostPtr[0], 1, dataLength);
    memset(&tmp_b.hostPtr[0], 2, dataLength);
    vecAddCPU(tmp_a.hostPtr, tmp_b.hostPtr, tmp_c.hostPtr, tmp_c.length);
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a.hostPtr[i] = 1;
    clReleaseMemObject(tmp_a.deviceBuffer);
    clReleaseMemObject(tmp_b.deviceBuffer);
    clReleaseMemObject(tmp_c.deviceBuffer);
    return time;
//    return 0;
}

double mallocVecBench(size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    auto *tmp_a = (float *) malloc(dataLength);
    auto *tmp_b = (float *) malloc(dataLength);
    auto *tmp_c = (float *) malloc(dataLength);
    time = getCurrentTime() - time;
    memset(&tmp_a[0], 1, dataLength);
    memset(&tmp_a[0], 2, dataLength);
    vecAddCPU(tmp_a, tmp_b, tmp_c, dataLength / sizeof(float));
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a[i] = 1;
    free(tmp_a);
    free(tmp_b);
    free(tmp_c);
    return time;
//    return 0;
}

double callocVecBench(size_t dataLength) {
    double time = 0;
//    time = getCurrentTime();
    auto *tmp_a = (float *) calloc(dataLength, sizeof(float));
    auto *tmp_b = (float *) calloc(dataLength, sizeof(float));
    auto *tmp_c = (float *) calloc(dataLength, sizeof(float));
    time = getCurrentTime() - time;
    memset(&tmp_a[0], 1, dataLength * sizeof(float));
    memset(&tmp_b[0], 2, dataLength * sizeof(float));
    vecAddCPU(tmp_a, tmp_b, tmp_c, dataLength);
    time = getCurrentTime() - time;
//    for (int i = 0; i < dataLength/ sizeof(float); ++i) tmp_a[i] = 1;
    free(tmp_a);
    free(tmp_b);
    free(tmp_c);
    return time;
    return 0;
}


int main() {
    setCurThreadAffinity(7);
    auto *cl = new cl_wrapper();
    const size_t LOOP_TIME = 100;
    size_t dataLength = 1024 * 1024 / 4;
    int pageSize = getpagesize();
    printf("Page size: %d \n", pageSize);

    for (int j = 0; j < 12; ++j) {
        flushed_printf("Memory alloc benchmark, size = %lf MB... \n", (double) dataLength / 1024 / 1024);

        double extra_time = 0;
        double t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += mallocBench(dataLength);
        }
        flushed_printf("\tmalloc runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

        extra_time = 0;
        t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += callocBench(dataLength / sizeof(float));
        }
        flushed_printf("\tcalloc runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

//        extra_time = 0;
//        t = getCurrentTime();
//        for (int i = 0; i < LOOP_TIME; ++i) {
//            extra_time += mmapBench(dataLength);
//        }
//        flushed_printf("\tmmap runtime: %lf sec ", getCurrentTime() - t);
//        flushed_printf("\t\tmmap extra_time = %lf s\n", extra_time);

        extra_time = 0;
        t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += zeroCopyBench(cl, dataLength);
        }
        flushed_printf("\tzeCopy runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

        dataLength = dataLength * 2;
    }

    dataLength = 1024 * 1024 / 4;
    for (int j = 0; j < 12; ++j) {
        flushed_printf("Memory alloc benchmark, size = %lf MB... \n", (double) dataLength / 1024 / 1024);

        double extra_time = 0;
        double t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += mallocVecBench(dataLength);
        }
        flushed_printf("\tmalloc vec runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

        extra_time = 0;
        t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += callocVecBench(dataLength / sizeof(float));
        }
        flushed_printf("\tcalloc vec runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

//        extra_time = 0;
//        t = getCurrentTime();
//        for (int i = 0; i < LOOP_TIME; ++i) {
//            extra_time += mmapBench(dataLength);
//        }
//        flushed_printf("\tmmap runtime: %lf sec ", getCurrentTime() - t);
//        flushed_printf("\t\tmmap extra_time = %lf s\n", extra_time);

        extra_time = 0;
        t = getCurrentTime();
        for (int i = 0; i < LOOP_TIME; ++i) {
            extra_time += zeroCopyVecBench(cl, dataLength);
        }
        flushed_printf("\tzeCopy vec runtime: %lf sec ", getCurrentTime() - t);
        flushed_printf("\textra_time = %lf s\n", extra_time);

        dataLength = dataLength * 2;
    }

    return 0;
}
