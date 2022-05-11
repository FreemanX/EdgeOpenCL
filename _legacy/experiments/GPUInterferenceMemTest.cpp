//
// Created by pfxu on 3/29/19.
//
#include "heteroCompLib.h"
#include <semaphore.h>
#include <pthread.h>
#include "timestamp.h"
#include <cstdio>
#include <alloca.h>

const char *KERNEL = "\n\
#define CMD_HELPER(FUNC, NAME) FUNC ## _ ## NAME\n\
#define CMD(FUNC, NAME) CMD_HELPER(FUNC, NAME)\n\
\n\
#define STRIDE (1 << STRIDE_ORDER)\n\
#define GRANULARITY (1 << GRANULARITY_ORDER) \n\
\n\
int reduce_int(int v)    { return v; }\n\
int reduce_int2(int2 v)  { return v.x+v.y; }\n\
int reduce_int4(int4 v)  { return v.x+v.y+v.z+v.w; }\n\
int reduce_int8(int8 v)  { return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7; }\n\
int reduce_int16(int16 v){ return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7+v.s8+v.s9+v.sA+v.sB+v.sC+v.sD+v.sE+v.sF; }\n\
\n\
__kernel void initialize(__global DATATYPE *data, const uint n, const int v) {\n\
	const uint id = get_global_id(0);\n\
	const uint stride = get_global_size(0);\n\
	uint offset = 0;\n\
	while( id+offset<n ){\n\
		data[id+offset] = (DATATYPE)v;\n\
		offset += stride;\n\
	}\n\
}\n\
\n\
__kernel void kernel1(__global DATATYPE *data, const uint n) {\n\
	// Get our global thread ID\n\
	const uint id = get_global_id(0);\n\
	const uint  low_order_id = id & ( (uint)(STRIDE - 1));\n\
	const uint high_order_id = id & (~(uint)(STRIDE - 1));\n\
	const uint index = (high_order_id << GRANULARITY_ORDER) | low_order_id;\n\
	const int localid = get_local_id(0);\n\
	const int group_size = get_local_size(0);\n\
/*	__local uint lcount;\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	if( localid==0 )\n\
		lcount = 0;*/\n\
	// Make sure we do not go out of bounds\n\
//	DATATYPE tmp = (DATATYPE)0;\n\
//	#pragma unroll\n\
	for(int i=0; i<GRANULARITY; i++){\n\
		DATATYPE v = data[index+i*STRIDE];\n\
		int v_sum = CMD(reduce, DATATYPE)(v);\n\
		if( v_sum>0 )\n\
			data[index+i*STRIDE] = (DATATYPE)0;\n\
		//tmp = tmp+data[index+i*STRIDE];\n\
		//data[index+i*STRIDE] = (DATATYPE)(index+i*STRIDE);\n\
	}\n\
/*	if( count )\n\
		atomic_add(&lcount, count);\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	// Atomic reduce in global memory\n\
	if(localid==0 && lcount){\n\
		atomic_add((__global int*)result, lcount);\n\
	}*/\n\
}\n\
";

void aligned_block_copy(int32_t *__restrict dst_,
                        int32_t *__restrict src,
                        long size) {
    volatile int32_t *dst = dst_;
    register int32_t t1, t2, t3, t4;
    register int32_t t5, t6, t7, t8;
    while ((size -= 64) >= 0) {
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        t5 = *src++;
        t6 = *src++;
        t7 = *src++;
        t8 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
        *dst++ = t5;
        *dst++ = t6;
        *dst++ = t7;
        *dst++ = t8;
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        t5 = *src++;
        t6 = *src++;
        t7 = *src++;
        t8 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
        *dst++ = t5;
        *dst++ = t6;
        *dst++ = t7;
        *dst++ = t8;
    }
}

void aligned_block_copy_backwards(int32_t *__restrict dst_,
                                  int32_t *__restrict src,
                                  long size) {
    volatile int32_t *dst = dst_;
    int32_t t1, t2, t3, t4;
    int32_t t5, t6, t7, t8;
    src += size / 4 - 1;
    dst += size / 4 - 1;
    while ((size -= 64) >= 0) {
        t1 = *src--;
        t2 = *src--;
        t3 = *src--;
        t4 = *src--;
        t5 = *src--;
        t6 = *src--;
        t7 = *src--;
        t8 = *src--;
        *dst-- = t1;
        *dst-- = t2;
        *dst-- = t3;
        *dst-- = t4;
        *dst-- = t5;
        *dst-- = t6;
        *dst-- = t7;
        *dst-- = t8;
        t1 = *src--;
        t2 = *src--;
        t3 = *src--;
        t4 = *src--;
        t5 = *src--;
        t6 = *src--;
        t7 = *src--;
        t8 = *src--;
        *dst-- = t1;
        *dst-- = t2;
        *dst-- = t3;
        *dst-- = t4;
        *dst-- = t5;
        *dst-- = t6;
        *dst-- = t7;
        *dst-- = t8;

    }
}

void random_read(int32_t *buffer, int32_t *no_use_submit_nullptr, long nbits) {
    uint32_t seed = 0;
    uintptr_t addrmask = (1 << nbits) - 1;
    uint32_t v;
    static volatile uint32_t dummy;
    //uint32_t rand1 = randNumberGeneratorCXX11(1, 123456);
    //uint32_t rand2 = randNumberGeneratorCXX11(1000000000, 2000000000);
    for (auto i = 0; i < 1000000; ++i) {
        seed = seed * 1103515245 + 12345;
        v = (seed >> 16) & 0xFF;
        seed = seed * 1103515245 + 12345;
        v |= (seed >> 8) & 0xFF00;
        seed = seed * 1103515245 + 12345;
        v |= seed & 0x7FFF0000;
        seed |= buffer[v & addrmask];

        //seed = seed * rand2 + rand1;
        //v = (seed >> 16) & 0xFF;
        //seed = seed * rand2 + rand1;
        //v |= (seed >> 8) & 0xFF00;
        //seed = seed * rand2 + rand1;
        //v |= seed & 0x7FFF0000;
        //seed |= buffer[v & addrmask];
    }

    dummy = seed;
}


void clmem_bench(
        unsigned int log2_indexes,
        unsigned int log2_grid,
        unsigned int log2_wgroup,
        unsigned int vecsize,
        unsigned int max_log2_stride,
        cl_kernel *kernels1,
        ZeroCopyMem<cl_int> *zeroCopyMem,
        cl_int index_space,
        cl_command_queue *exeQueue) {
    const size_t glWS[1] = {index_space / pow2(log2_indexes - log2_grid)};
    const size_t lcWS[1] = {pow2(log2_wgroup)};
    auto *total_times = (double *) alloca(sizeof(double) * (max_log2_stride + 1));
    const int REPETITIONS = 16;
    flushed_printf("Running... \n");
    cl_event ev_wait;
    for (int stride_offset = 0; stride_offset <= (int) max_log2_stride; stride_offset++) {
        clSetKernelArg(kernels1[stride_offset], 0, sizeof(cl_mem), &zeroCopyMem->deviceBuffer);
        clSetKernelArg(kernels1[stride_offset], 1, sizeof(cl_int), &index_space);
        clEnqueueNDRangeKernel(*exeQueue, kernels1[stride_offset], 1, nullptr, glWS, lcWS, 0, nullptr, nullptr);
        clFinish(*exeQueue);
        // run benchmarks multiple times
        auto *times = (double *) alloca(sizeof(double) * REPETITIONS);
        clSetKernelArg(kernels1[stride_offset], 0, sizeof(cl_mem), &zeroCopyMem->deviceBuffer);
        clSetKernelArg(kernels1[stride_offset], 1, sizeof(cl_int), &index_space);
        int do_run_experiment = 1;
        while (do_run_experiment) {
            for (int i = 0; i < REPETITIONS; i++) {
                clEnqueueNDRangeKernel(*exeQueue, kernels1[stride_offset], 1, nullptr, glWS, lcWS, 0, nullptr,
                                       &ev_wait);
                clWaitForEvents(1, &ev_wait);
                times[i] = cl_wrapper::getExecutionTime(&ev_wait);
            }
            qsort(times, REPETITIONS, sizeof(times[0]), compare_doubles);
            const double median_time = REPETITIONS % 2 ? times[REPETITIONS / 2] :
                                       (times[REPETITIONS / 2 - 1] + times[REPETITIONS / 2]) / 2;
            double average_time = 0., variance = 0.;
            for (int i = 0; i < REPETITIONS; i++)
                average_time += times[i];
            average_time /= REPETITIONS;
            for (int i = 0; i < REPETITIONS; i++)
                variance += sqr(times[i] - average_time);
            variance /= REPETITIONS;
            double variation_coeff = sqrt(variance) / average_time;
            const double VAR_COEFF_THRESHOLD = 0.3;
            if (variation_coeff < VAR_COEFF_THRESHOLD)
                do_run_experiment = 0;
            total_times[stride_offset] = median_time;//average_time;
        }
    }

    printf("\nSummary:");
    for (int stride_offset = 0; stride_offset <= (int) max_log2_stride; stride_offset++) {
        printf("\nStride magnitude %2d: Bandwidth %7.3f GB/sec (avg time %10f msecs)",
               stride_offset,
               pow2(log2_indexes) * vecsize * sizeof(int) / (total_times[stride_offset] * 1000.0 * 1000.0 * 1000.0),
               1000.0 * total_times[stride_offset]);
        if (stride_offset == log2_grid) printf(" *Grid striding");
        if (stride_offset == log2_wgroup) printf(" *Workgroup striding");
        if (stride_offset == 0 && max_log2_stride > 0) printf(" *Serial accesses");
    }
    printf("\n");
}

typedef struct {
    int32_t *dst;
    int32_t *src;
    long size; // in byte
    int cpu_id;

    void (*f)(int32_t *, int32_t *, long);

    sem_t *semSync;
    sem_t *mainSyncSem;
} thread_data;

void *interference_memory(void *arg) {
    auto *threadData = (thread_data *) arg;
    //if (setCurThreadAffinity(threadData->cpu_id)) {
    //    checkCurThreadAffinity(0);
    //    flushed_printf("pthread launch on core %d failed!\n", threadData->cpu_id);
    //    pthread_exit(nullptr);
    //}
    while (setCurThreadAffinity(threadData->cpu_id) != 0) {//
        //  try to set affinity, a war must win!
    }
    if(threadData->mainSyncSem) sem_post(threadData->mainSyncSem);
    flushed_printf("Memory interference thread... %d\n", threadData->cpu_id);

    while (true) {
        int32_t *dst_ = threadData->dst;
        int32_t *src_ = threadData->src;
        long size_ = threadData->size;
        threadData->f(dst_, src_, size_);
        if (sem_trywait(threadData->semSync) == 0) break;
    }

    flushed_printf("Memory interference done...%d\n", threadData->cpu_id);
    pthread_exit(nullptr);
}

int main(int argc, char **argv) {
    //Bench Config
    unsigned int log2_indexes = 26;
    unsigned int log2_grid = 20;
    unsigned int log2_wgroup = 8;
    unsigned int vecsize = 2;
    const unsigned int max_log2_stride = log2_indexes > log2_grid ? log2_grid : 0;
    char s_vecsize[3] = "";
    if (vecsize > 1)
        sprintf(s_vecsize, "%d", vecsize);

    printf("\nBenchmark parameters:\n");
    printf("log2_indexes = %d, log2_grid = %d, log2_wgroup = %d, vecsize = %d \n",
           log2_indexes, log2_grid, log2_wgroup, vecsize);
    printf("index space     : %d\n", pow2(log2_indexes));
    printf("vector width    : %d (type: int%s)\n", vecsize, s_vecsize);
    //printf("element space   : %d\n", pow2(log2_indexes)*vecsize);
    {
        unsigned long int req_mem = (unsigned long int) (pow2(log2_indexes)) * vecsize * sizeof(int) / 1024;
        char req_mem_unit = 'K';
        if (req_mem >= 1024 * 8) {
            req_mem /= 1024;
            req_mem_unit = 'M';
        }
        if (req_mem >= 1024 * 8) {
            req_mem /= 1024;
            req_mem_unit = 'G';
        }
        printf("required memory : %lu %cB\n", req_mem, req_mem_unit);
    }
    printf("grid space      : %d (%d workgroups)\n", pow2(log2_grid), pow2(log2_grid - log2_wgroup));
    printf("workgroup size  : %d\n", pow2(log2_wgroup));
    //printf("total workgroups: %d\n", pow2(log2_grid-log2_wgroup));
    printf("granularity     : %d\n", pow2(log2_indexes - log2_grid));

    setCurThreadAffinity(0);
    checkCurThreadAffinity(0);

    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    auto *cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    char options[4096];
    auto *programs = (cl_program *) alloca(sizeof(cl_program) * (max_log2_stride + 1));
    for (int stride_offset = 0; stride_offset <= (int) max_log2_stride; stride_offset++) {
        sprintf(options, " -DDATATYPE=int%s -DSTRIDE_ORDER=%d -DGRANULARITY_ORDER=%d",
                s_vecsize, stride_offset, log2_indexes - log2_grid);
        programs[stride_offset] = cl->createProgramWithOptions(&KERNEL, 1, options);
    }
    auto *kernels_init = (cl_kernel *) alloca(sizeof(cl_kernel) * (max_log2_stride + 1));
    auto *kernels1 = (cl_kernel *) alloca(sizeof(cl_kernel) * (max_log2_stride + 1));
    for (int stride_offset = 0; stride_offset <= (int) max_log2_stride; stride_offset++) {
        kernels_init[stride_offset] = cl->createKernel("initialize", programs[stride_offset]);
        kernels1[stride_offset] = cl->createKernel("kernel1", programs[stride_offset]);
    }
    cl_command_queue exeQueue = cl->createProfilingCmdQueue();

    cl_int index_space = pow2(log2_indexes);
    const cl_int zero = 0;

    flushed_printf("Requested buffer size: %d MB\n", index_space * vecsize * sizeof(cl_int) / 1024 / 1024);
    auto zeroCopyMem = init_zero_copy_region<cl_int>(cl, index_space * vecsize);
    const size_t glWS[1] = {index_space / pow2(log2_indexes - log2_grid)};
    const size_t lcWS[1] = {pow2(log2_wgroup)};

    flushed_printf("Zeroing buffer... \n");
    clSetKernelArg(kernels_init[0], 0, sizeof(cl_mem), &zeroCopyMem.deviceBuffer);
    clSetKernelArg(kernels_init[0], 1, sizeof(cl_int), &index_space);
    clSetKernelArg(kernels_init[0], 2, sizeof(cl_int), &zero);
    clEnqueueNDRangeKernel(exeQueue, kernels_init[0], 1, nullptr, glWS, lcWS, 0, nullptr, nullptr);

    sem_t syncSem;
    sem_init(&syncSem, 0, 0);

    flushed_printf("\n\n======No interference====\n");
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);

    thread_data threadData = {
            .dst     = zeroCopyMem.hostPtr,
            .src     = zeroCopyMem.hostPtr + (index_space * vecsize / 2) - 1,
            .size    = static_cast<long>(index_space * vecsize * sizeof(cl_int) / 2),
            .cpu_id  = 0,
            .f       = aligned_block_copy,
            .semSync = &syncSem};
    int rc;
    flushed_printf("\n\n======Forward interference Little====\n");
    pthread_t fthread_l;
    if ((rc = pthread_create(&fthread_l, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(fthread_l, nullptr);

    flushed_printf("\n\n======Forward interference Big====\n");
    pthread_t fthread_b;
    threadData.cpu_id = 7;
    if ((rc = pthread_create(&fthread_b, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(fthread_b, nullptr);

    flushed_printf("\n\n======Multiple CPU cores interference ====\n");
    pthread_t mc_threads[8];
    thread_data data[8];
    sem_t syncSems[8];
    sem_t mainSyncSems[8];

    for (int i = 0; i < 8; ++i) {
        sem_init(&syncSems[i],     0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        data[i] = {
                .dst     = zeroCopyMem.hostPtr,
                .src     = zeroCopyMem.hostPtr + (index_space * vecsize / 2) - 1,
                .size    = static_cast<long>(index_space * vecsize * sizeof(cl_int) / 2),
                .cpu_id  = i,
                .f       = aligned_block_copy,
                .semSync = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&mc_threads[i], nullptr, interference_memory, &data[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 0; i < 8; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    for (int i = 0; i < 8; i++) {
        sem_post(&syncSems[i]);
        pthread_join(mc_threads[i], nullptr);
    }

    flushed_printf("\n\n======Multiple CPU cores random access interference ====\n");
    for (int i = 0; i < 8; ++i) {
        sem_init(&syncSems[i],     0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        data[i] = {
                .dst     = zeroCopyMem.hostPtr,
                .src     = zeroCopyMem.hostPtr + (index_space * vecsize / 2) - 1,
                .size    = static_cast<long>(index_space * vecsize * sizeof(cl_int) / 2),
                .cpu_id  = i,
                .f       = random_read,
                .semSync = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&mc_threads[i], nullptr, interference_memory, &data[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 0; i < 8; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    for (int i = 0; i < 8; i++) {
        sem_post(&syncSems[i]);
        pthread_join(mc_threads[i], nullptr);
    }


    flushed_printf("\n\n======Backward interference Little====\n");
    pthread_t bthread_l;
    threadData.cpu_id = 0;
    threadData.f = aligned_block_copy_backwards;
    if ((rc = pthread_create(&bthread_l, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(bthread_l, nullptr);

    flushed_printf("\n\n======Backward interference Big====\n");
    pthread_t bthread_b;
    threadData.cpu_id = 7;
    if ((rc = pthread_create(&bthread_b, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(bthread_b, nullptr);

    flushed_printf("\n\n======Random Read interference Little====\n");
    pthread_t rthread_l;
    threadData.cpu_id = 0;
    threadData.f = random_read;
    threadData.dst = zeroCopyMem.hostPtr;
    threadData.src = nullptr;
    threadData.size = log2_indexes + (long) log2(vecsize);
    if ((rc = pthread_create(&rthread_l, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(rthread_l, nullptr);

    flushed_printf("\n\n======Random Read interference Big====\n");
    pthread_t rthread_b;
    threadData.cpu_id = 7;
    if ((rc = pthread_create(&rthread_b, nullptr, interference_memory, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    clmem_bench(log2_indexes, log2_grid, log2_wgroup, vecsize,
                max_log2_stride, kernels1, &zeroCopyMem, index_space, &exeQueue);
    sem_post(&syncSem);
    pthread_join(rthread_b, nullptr);

    sem_destroy(&syncSem);
    for (int i = 0; i < 8; ++i) {
        sem_destroy(&syncSems[i]);
        sem_destroy(&mainSyncSems[i]);
    }
    clReleaseMemObject(zeroCopyMem.deviceBuffer);
    return 0;
}