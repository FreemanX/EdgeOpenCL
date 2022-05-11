//
// Created by pfxu on 5/24/19.
//

#include "heteroCompLib.h"
#include <algorithm>

using namespace std;

const char *kernelSource = "\n" \
"__kernel void vecAdd(  __global float *a,                       \n" \
"                       __global float *b,                       \n" \
"                       __global float *c,                       \n" \
"                       const int n)                             \n" \
"{                                                               \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    if (id < n)                                                 \n" \
"        c[id] = a[id] + b[id];                                  \n" \
"}                                                               \n";

cl_wrapper *cl;
cl_program program;
cl_kernel kernel;
cl_command_queue exeCmdQueues[8];


struct Task_pack {
    size_t offset;
    size_t workingLength;
};

struct Worker_pack {
    int proc_id;
    ZeroCopyMem<float> *a;
    ZeroCopyMem<float> *b;
    ZeroCopyMem<float> *c;
    ConcurrentQ<Task_pack> **taskQueues;
    int numQ;
    float *work_cnt;
};

void initCLEnv() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices->deviceId);
    program = cl->createProgram(&kernelSource, 1);
    kernel = cl->createKernel("vecAdd", program);
    for (auto &exeCmdQueue : exeCmdQueues) {
        exeCmdQueue = cl->createProfilingCmdQueue();
    }
}

void CPUMatAdd(float *a, float *b, float *c, int dataLength) {
    for (int i = 0; i < dataLength; ++i) {
        c[i] = a[i] + b[i];
    }
}

void printTaskPackQueues(ConcurrentQ<Task_pack> **queues, int numQ, size_t dataLength) {
    std::cout << "======Data length: " << dataLength << ", task queues======\n";
    for (int q = 0; q < numQ; ++q) {
        std::cout << "==== task queue" << q << "====\n";
        std::cout << "idx\toffset\tlength\to+t\n";
        for (int i = 0; i < queues[q]->size; ++i) {
            std::cout << i << "\t" << queues[q]->array[(queues[q]->front + i) % queues[q]->capacity].offset << "\t";
            std::cout << queues[q]->array[(queues[q]->front + i) % queues[q]->capacity].workingLength << "\t";
            std::cout << queues[q]->array[(queues[q]->front + i) % queues[q]->capacity].workingLength +
                         queues[q]->array[(queues[q]->front + i) % queues[q]->capacity].offset << "\n";
        }
        std::cout << "(size " << queues[q]->size;
        std::cout << ", front " << queues[q]->front;
        std::cout << ", rear " << queues[q]->rear;
        std::cout << ")\n";
    }
    std::cout << "\n\n";
}

bool queuesAllEmpty(ConcurrentQ<Task_pack> **queues, int numQ) {
    bool allEmpty = true;
    for (int i = 0; i < numQ; ++i) {
        allEmpty &= ConcurrentQ_isEmpty(queues[i]);
    }
    return allEmpty;
}

void *Worker(void *worker_pack) {
    auto *workerPack = (Worker_pack *) worker_pack;
    Task_pack taskPack{};
    ConcurrentQ<Task_pack> *taskQueue;

    if (workerPack->proc_id < 8) { //CPU
        while (setCurThreadAffinity(workerPack->proc_id) != 0) {
            if (queuesAllEmpty(workerPack->taskQueues, workerPack->numQ)) pthread_exit(nullptr);
        }
        for (int i = 0; i < workerPack->numQ; ++i) {
            taskQueue = workerPack->taskQueues[(workerPack->proc_id + i) % workerPack->numQ];
            while (!ConcurrentQ_isEmpty(taskQueue)) {
                if (ConcurrentQ_dequeueF<Task_pack>(taskQueue, &taskPack)) {
                    CPUMatAdd(workerPack->a->hostPtr + taskPack.offset,
                              workerPack->b->hostPtr + taskPack.offset,
                              workerPack->c->hostPtr + taskPack.offset,
                              taskPack.workingLength);
                    workerPack->work_cnt[workerPack->proc_id]++;
                }
            }
        }
    } else { //GPU
        setCurThreadAffinity((workerPack->proc_id + 1) % 4);
        cl_wrapper::checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), &workerPack->a->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 1, sizeof(cl_mem), &workerPack->b->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 2, sizeof(cl_mem), &workerPack->c->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 3, sizeof(int), (int *) &workerPack->a->length));
        cl_event e[8];
        for (int i = 0; i < workerPack->numQ; ++i) {
            taskQueue = workerPack->taskQueues[(workerPack->proc_id + i) % workerPack->numQ];
            while (!ConcurrentQ_isEmpty(taskQueue)) {
                int eventCnt = 0;
                for (int j = 0; j < 4; ++j) { // TODO: have to lock 4 times, how about lock once and get 4 packs
                    if (ConcurrentQ_dequeueF<Task_pack>(taskQueue, &taskPack)) {
                        size_t global_size[] = {taskPack.workingLength};
                        size_t offset[] = {taskPack.offset};
                        cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueues[j],
                                                                      kernel,
                                                                      1, offset,
                                                                      global_size,
                                                                      nullptr, 0, nullptr, &e[j]));
                        workerPack->work_cnt[workerPack->proc_id]++;
                        eventCnt = i + 1;
                    }
                }

                for (int j = 0; j < eventCnt; ++j) {
                    clWaitForEvents(1, &e[j]);
                }

            }
        }
    }
    pthread_exit(nullptr);
}

void executionLayer(ZeroCopyMem<float> *vecA,
                    ZeroCopyMem<float> *vecB,
                    ZeroCopyMem<float> *vecC,
                    float *workCnt,
                    int s,
                    int e,
                    ConcurrentQ<Task_pack> **taskQueues,
                    int numQ) {
    int numWorker = e - s;
    pthread_t workers[numWorker];
    Worker_pack workerPacks[numWorker];

    for (int i = 0; i < numWorker; ++i) {
        workerPacks[i].a = vecA;
        workerPacks[i].b = vecB;
        workerPacks[i].c = vecC;
        workerPacks[i].proc_id = i + s;
        workerPacks[i].taskQueues = taskQueues;
        workerPacks[i].work_cnt = workCnt;
        workerPacks[i].numQ = numQ;
        if (pthread_create(&workers[i], nullptr, &Worker, &workerPacks[i])) {
            flushed_printf("Pthread %d create failed...\n", workerPacks[i].proc_id);
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        }
    }

    for (int i = 0; i < numWorker; ++i) {
        pthread_join(workers[i], nullptr);
    }
}

void genTaskQueue(int numQ, int numParts, int offsetIdx, size_t dataLength,
                  ConcurrentQ<Task_pack> **queues) {
    Task_pack taskPack{};
    size_t workingLength = dataLength / numParts;
    taskPack.workingLength = workingLength;
    unsigned capacity = std::max(numParts / numQ, 1);
    for (int q = 0; q < numQ; ++q) {
        queues[q] = ConcurrentQ_create<Task_pack>(capacity + 1);
        for (int i = 0; i <= capacity; ++i) {
            taskPack.offset = offsetIdx;
            if (dataLength - offsetIdx < 2 * workingLength) {
                if ((taskPack.workingLength = dataLength - taskPack.offset) > 0) {
                    ConcurrentQ_enqueue(queues[q], taskPack);
                    offsetIdx += workingLength;
                }
                break;
            }
            ConcurrentQ_enqueue(queues[q], taskPack);
            offsetIdx += workingLength;
        }
    }
}

void genTaskQueuesEven(int numQ, int numParts, size_t dataLength,
                       ConcurrentQ<Task_pack> **queues) {
    genTaskQueue(numQ, numParts, 0, dataLength, queues);
}

void genTaskQueuesUneven(int numQ, float GPUSplit, int numParts, size_t dataLength,
                         ConcurrentQ<Task_pack> **queues) {
    if (GPUSplit > 0) {
        //TODO: numQ = 1; numParts < numQ
        size_t cpuEnd = std::floor(dataLength * (1 - GPUSplit));
//        genTaskQueue(numQ - 1, numParts, 0, cpuEnd, queues);
//        genTaskQueue(1, 1, cpuEnd, dataLength, queues + numQ - 1);
        size_t cpuNumParts = std::floor(numParts * (1 - GPUSplit));
        genTaskQueue(numQ - 1, cpuNumParts, 0, cpuEnd, queues);
        genTaskQueue(1, numParts - cpuNumParts, cpuEnd, dataLength, queues + numQ - 1);
    } else {
        genTaskQueuesEven(numQ, numParts, dataLength, queues);
    }
}


void scheduleLayer(ZeroCopyMem<float> *a,
                   ZeroCopyMem<float> *b,
                   ZeroCopyMem<float> *c) {
    auto *c_cpu = (float *) malloc(a->length * sizeof(float));
    vecAddCPU(a->hostPtr, b->hostPtr, c_cpu, a->length);

    int LOOP_TIME = 10;
    float *workCnt = (float *) malloc(20 * sizeof(float));
    int s, e;
    int numParts;
    double timer;

    s = 8, e = 9;
    cl_event clEvent;
    cl_int cl_zero = 0;
    cl_command_queue exeCmdQueue = cl->createProfilingCmdQueue();
    cout << "Processors s: " << s << ", e: " << e << "\n";
    int numQ = e - s;
    timer = getCurrentTime();
    for (int j = 0; j < 14; ++j) {
        memset(workCnt, 0, 20 * sizeof(float));
        numParts = pow2(j);
        auto **queues = (ConcurrentQ<Task_pack> **) malloc(numQ * sizeof(ConcurrentQ<Task_pack> *));
        for (int i = 0; i < LOOP_TIME; ++i) {
            genTaskQueuesEven(numQ, numParts, a->length, queues);
            executionLayer(a, b, c, workCnt, s, e, queues, numQ);
        }
//        for (int i = 0; i < LOOP_TIME; ++i)


        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
        printf("\n");

        cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
        clWaitForEvents(1, &clEvent);
        free(queues);
    }
}

int main(int argc, char **argv) {

    const size_t w = 1024;
    const size_t h = 1024 / 4;
    const size_t dataLength = w * h * 512 / sizeof(float);
    initCLEnv();
    double timer = getCurrentTime();
    ZeroCopyMem<float> a = init_zero_copy_region<float>(cl, dataLength);
    ZeroCopyMem<float> b = init_zero_copy_region<float>(cl, dataLength);
    randArrayGenerator<float>(0.0f, 1.0, a.hostPtr, dataLength);
    randArrayGenerator<float>(0.0f, 1.0, b.hostPtr, dataLength);
    ZeroCopyMem<float> c = init_zero_copy_region<float>(cl, dataLength);
    cout << "Init time: " << getCurrentTime() - timer << endl;
    timer = getCurrentTime();
    scheduleLayer(&a, &b, &c);
    timer = getCurrentTime() - timer;

    cout << "Total time: " << timer << "s \n";


//    int numQ = 5;
//    int numParts = 4;
//    size_t dataLength = 1024 * 1000;
//    auto **queues = (ConcurrentQ<Task_pack> **) malloc(numQ * sizeof(ConcurrentQ<Task_pack> *));
////    genTaskQueuesUneven(numQ, 0.6, numParts, dataLength, queues);
//    genTaskQueuesEven(numQ, numParts, dataLength, queues);
//    printTaskPackQueues(queues, numQ, dataLength);
//
//    free(queues);
//    numQ = 1;
//    queues = (ConcurrentQ<Task_pack> **) malloc(numQ * sizeof(ConcurrentQ<Task_pack> *));
////    genTaskQueuesUneven(numQ, 0.6, numParts, dataLength, queues);
//    genTaskQueuesEven(numQ, numParts, dataLength, queues);
//    printTaskPackQueues(queues, numQ, dataLength);
//
//    free(queues);
//    numQ = 2;
//    queues = (ConcurrentQ<Task_pack> **) malloc(numQ * sizeof(ConcurrentQ<Task_pack> *));
////    genTaskQueuesUneven(numQ, 0.6, numParts, dataLength, queues);
//    genTaskQueuesEven(numQ, numParts, dataLength, queues);
//    printTaskPackQueues(queues, numQ, dataLength);
//
//    free(queues);
//    numQ = 3;
//    queues = (ConcurrentQ<Task_pack> **) malloc(numQ * sizeof(ConcurrentQ<Task_pack> *));
////    genTaskQueuesUneven(numQ, 0.6, numParts, dataLength, queues);
//    genTaskQueuesEven(numQ, numParts, dataLength, queues);
//    printTaskPackQueues(queues, numQ, dataLength);
    return 0;
}
