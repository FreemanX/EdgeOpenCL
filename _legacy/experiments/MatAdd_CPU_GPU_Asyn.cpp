//
// Created by pfxu on 5/15/19.
//
//Prototyping hetero-compute concurrent mode


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
"        c[id] = a[id] + a[id] * b[id];                                  \n" \
"}                                                               \n";

cl_wrapper *cl;
cl_program program;
cl_kernel kernel;
cl_command_queue exeCmdQueues[8];

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
        c[i] = a[i] + a[i] * b[i];
    }
}

struct Task_pack {
    size_t offset;
    size_t workingLength;
};

struct Worker_pack {
    int proc_id;
    ZeroCopyMem<float> *a;
    ZeroCopyMem<float> *b;
    ZeroCopyMem<float> *c;
    ConcurrentQ<Task_pack> *taskQueue;
    float *work_cnt;
};

void *HardWorker(void *worker_pack) {
    auto *workerPack = (Worker_pack *) worker_pack;
    int numPack = 1;
    auto *taskPacks = (Task_pack *) malloc(numPack * sizeof(Task_pack));
    int n;
    if (workerPack->proc_id < 8) { //CPU
        while (setCurThreadAffinity(workerPack->proc_id) != 0) {
            if (ConcurrentQ_isEmpty(workerPack->taskQueue)) pthread_exit(nullptr);
        }
        while (!ConcurrentQ_isEmpty(workerPack->taskQueue)) {
            if ((n = ConcurrentQ_dequeueN<Task_pack>(workerPack->taskQueue, taskPacks, numPack)) > 0) {
                CPUMatAdd(workerPack->a->hostPtr + taskPacks[0].offset,
                          workerPack->b->hostPtr + taskPacks[0].offset,
                          workerPack->c->hostPtr + taskPacks[0].offset,
                          taskPacks[0].workingLength * n);
                workerPack->work_cnt[workerPack->proc_id] += n;
            }
            numPack++;
            auto tmp_ptr = static_cast<Task_pack *>(realloc(taskPacks, numPack * sizeof(Task_pack)));
            if (tmp_ptr) taskPacks = tmp_ptr;
            else numPack--;
        }
    } else { //GPU
        setCurThreadAffinity((workerPack->proc_id + 1) % 4);
        cl_wrapper::checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), &workerPack->a->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 1, sizeof(cl_mem), &workerPack->b->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 2, sizeof(cl_mem), &workerPack->c->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 3, sizeof(int), (int *) &workerPack->a->length));
        cl_event e[8];
        int eventCnt = 0;
        while (!ConcurrentQ_isEmpty(workerPack->taskQueue)) {
            for (int i = 0; i < 4; ++i) { // TODO: have to lock 4 times, how about lock once and get 4 packs
                if ((n = ConcurrentQ_dequeueN<Task_pack>(workerPack->taskQueue, taskPacks, numPack)) > 0) {
                    size_t global_size[] = {taskPacks[0].workingLength * n};
                    size_t offset[] = {taskPacks[0].offset};
                    cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueues[i],
                                                                  kernel,
                                                                  1, offset,
                                                                  global_size,
                                                                  nullptr, 0, nullptr, &e[i]));
//                                                              wg_size, 0, nullptr, nullptr));
                    workerPack->work_cnt[workerPack->proc_id] += n;
                    eventCnt = i + 1;
                }
            }
            numPack *= 8;
            auto tmp_ptr = static_cast<Task_pack *>(realloc(taskPacks, numPack * sizeof(Task_pack)));
            if (tmp_ptr) taskPacks = tmp_ptr;
            else numPack /= 8;
            for (int j = 0; j < eventCnt; ++j) {
                clWaitForEvents(1, &e[j]);
            }
        }
//        clReleaseEvent(e);
//        std::cout << "\n";
    }

    pthread_exit(nullptr);
}


void *Worker(void *worker_pack) {
    auto *workerPack = (Worker_pack *) worker_pack;
    Task_pack taskPack{};
    if (workerPack->proc_id < 8) { //CPU
        while (setCurThreadAffinity(workerPack->proc_id) != 0) {
            if (ConcurrentQ_isEmpty(workerPack->taskQueue)) pthread_exit(nullptr);
        }
        while (!ConcurrentQ_isEmpty(workerPack->taskQueue)) {
            if (ConcurrentQ_dequeueF<Task_pack>(workerPack->taskQueue, &taskPack)) {
                CPUMatAdd(workerPack->a->hostPtr + taskPack.offset,
                          workerPack->b->hostPtr + taskPack.offset,
                          workerPack->c->hostPtr + taskPack.offset,
                          taskPack.workingLength);
//                flushed_printf("core %d take pack %d\n", workerPack->proc_id, taskPack.offset);
                workerPack->work_cnt[workerPack->proc_id]++;
            }
        }
    } else { //GPU
//        setCurThreadAffinity(7);
        setCurThreadAffinity((workerPack->proc_id + 1) % 4);
        cl_wrapper::checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), &workerPack->a->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 1, sizeof(cl_mem), &workerPack->b->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 2, sizeof(cl_mem), &workerPack->c->deviceBuffer));
        cl_wrapper::checkError(clSetKernelArg(kernel, 3, sizeof(int), (int *) &workerPack->a->length));
//        cl_command_queue exeCmdQueue;
//        exeCmdQueue = cl->createProfilingCmdQueue();
        cl_event e[8];
        int eventCnt = 0;
        while (!ConcurrentQ_isEmpty(workerPack->taskQueue)) {
            for (int i = 0; i < 1; ++i) { // TODO: have to lock 4 times, how about lock once and get 4 packs
                if (ConcurrentQ_dequeueF<Task_pack>(workerPack->taskQueue, &taskPack)) {
                    size_t global_size[] = {taskPack.workingLength};
//                    size_t wg_size[] = {1024};
                    size_t offset[] = {taskPack.offset};
                    cl_wrapper::checkError(clEnqueueNDRangeKernel(exeCmdQueues[i],
                                                                  kernel,
                                                                  1, offset,
                                                                  global_size,
                                                                  nullptr, 0, nullptr, &e[i]));
//                                                              wg_size, 0, nullptr, nullptr));
                    workerPack->work_cnt[workerPack->proc_id]++;
                    eventCnt = i + 1;
                }
            }
            for (int j = 0; j < eventCnt; ++j) {
                clWaitForEvents(1, &e[j]);
            }
        }
//        clReleaseEvent(e);
//        std::cout << "\n";
    }

    pthread_exit(nullptr);
}

void executionLayer(ZeroCopyMem<float> *vecA,
                    ZeroCopyMem<float> *vecB,
                    ZeroCopyMem<float> *vecC,
                    float *workCnt,
                    int s,
                    int e,
                    ConcurrentQ<Task_pack> *taskQueue) {
    int numWorker = e - s;
    pthread_t workers[numWorker];
    Worker_pack workerPacks[numWorker];

    for (int i = 0; i < numWorker; ++i) {
        workerPacks[i].a = vecA;
        workerPacks[i].b = vecB;
        workerPacks[i].c = vecC;
        workerPacks[i].proc_id = i + s;
        workerPacks[i].taskQueue = taskQueue;
        workerPacks[i].work_cnt = workCnt;
//        flushed_printf("Execute on %d \n", workerPacks[i].proc_id);
//        if (pthread_create(&workers[i], nullptr, &HardWorker, &workerPacks[i])) {
        if (pthread_create(&workers[i], nullptr, &Worker, &workerPacks[i])) {
            flushed_printf("Pthread %d create failed...\n", workerPacks[i].proc_id);
            flushed_printf("Error %d: %s\n", errno, strerror(errno));
        }
    }

    for (int i = 0; i < numWorker; ++i) {
        pthread_join(workers[i], nullptr);
    }
}


void scheduleLayer(ZeroCopyMem<float> *a,
                   ZeroCopyMem<float> *b,
                   ZeroCopyMem<float> *c) {
    auto *c_cpu = (float *) malloc(a->length * sizeof(float));
    CPUMatAdd(a->hostPtr, b->hostPtr, c_cpu, a->length);

    unsigned numParts;
    int LOOP_TIME = 100;
    Task_pack taskPack{};
    float *workCnt = (float *) malloc(20 * sizeof(float));
    double timer = 0;
    int s, e;

    s = 8, e = 9;
    cl_event clEvent;
    cl_int cl_zero = 0;
    cl_command_queue exeCmdQueue = cl->createProfilingCmdQueue();
    cout << "Processors s: " << s << ", e: " << e << "\n";
    printf("%-12s%-10s\t%s\n", "numParts", "Time", "Work count");
    for (int j = 0; j < 14; ++j) {
        memset(workCnt, 0, 20 * sizeof(float));
        numParts = pow2(j);
        ConcurrentQ<Task_pack> *queue = ConcurrentQ_create<Task_pack>(numParts);
        timer = getCurrentTime();
        for (int k = 0; k < LOOP_TIME; ++k) {
            size_t offsetIdx = 0;
            for (int i = 0; i < numParts; ++i) {
                taskPack.workingLength = a->length / numParts;
                taskPack.offset = offsetIdx;
                offsetIdx += taskPack.workingLength;
                if (a->length - offsetIdx < taskPack.workingLength || i == numParts - 1) {
                    taskPack.workingLength = a->length - taskPack.offset;
                    ConcurrentQ_enqueue(queue, taskPack);
                    break;
                }
                ConcurrentQ_enqueue(queue, taskPack);
            }
            executionLayer(a, b, c, workCnt, s, e, queue);
        }
        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
        printf("\n");

        cl_wrapper::checkError(clEnqueueFillBuffer(cl->memCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
        clWaitForEvents(1, &clEvent);
        ConcurrentQ_destroy(queue);
    }

    s = 4, e = 9;
    cout << "Processors s: " << s << ", e: " << e << "\n";
    printf("%-12s%-10s\t%s\n", "numParts", "Time", "Work count");
    for (int j = 0; j < 14; ++j) {
        memset(workCnt, 0, 20 * sizeof(float));
        numParts = pow2(j);
        ConcurrentQ<Task_pack> *queue = ConcurrentQ_create<Task_pack>(numParts);
        timer = getCurrentTime();
        for (int k = 0; k < LOOP_TIME; ++k) {
            size_t offsetIdx = 0;
            for (int i = 0; i < numParts; ++i) {
                taskPack.workingLength = a->length / numParts;
                taskPack.offset = offsetIdx;
                offsetIdx += taskPack.workingLength;
                if (a->length - offsetIdx < taskPack.workingLength || i == numParts - 1) {
                    taskPack.workingLength = a->length - taskPack.offset;
                    ConcurrentQ_enqueue(queue, taskPack);
                    break;
                }
                ConcurrentQ_enqueue(queue, taskPack);
            }
            executionLayer(a, b, c, workCnt, s, e, queue);
        }
        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
        printf("\n");

        cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
        clWaitForEvents(1, &clEvent);
        ConcurrentQ_destroy(queue);
    }

    s = 4, e = 8;
    cout << "Processors s: " << s << ", e: " << e << "\n";
    printf("%-12s%-10s\t%s\n", "numParts", "Time", "Work count");
    for (int j = 0; j < 14; ++j) {
        memset(workCnt, 0, 20 * sizeof(float));
        numParts = pow2(j);
        ConcurrentQ<Task_pack> *queue = ConcurrentQ_create<Task_pack>(numParts);
        timer = getCurrentTime();
        for (int k = 0; k < LOOP_TIME; ++k) {
            size_t offsetIdx = 0;
            for (int i = 0; i < numParts; ++i) {
                taskPack.workingLength = a->length / numParts;
                taskPack.offset = offsetIdx;
                offsetIdx += taskPack.workingLength;
                if (a->length - offsetIdx < taskPack.workingLength || i == numParts - 1) {
                    taskPack.workingLength = a->length - taskPack.offset;
                    ConcurrentQ_enqueue(queue, taskPack);
                    break;
                }
                ConcurrentQ_enqueue(queue, taskPack);
            }
            executionLayer(a, b, c, workCnt, s, e, queue);
        }
        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
        printf("\n");

        cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
        clWaitForEvents(1, &clEvent);
        ConcurrentQ_destroy(queue);
    }

//    s = 0, e = 8;
//    cout << "Processors s: " << s << ", e: " << e << "\n";
//    printf("%-12s%-10s\t%s\n", "numParts", "Time", "Work count");
//    for (int j = 0; j < 14; ++j) {
//        memset(workCnt, 0, 20 * sizeof(float));
//        numParts = pow2(j);
//        ConcurrentQ<Task_pack> *queue = ConcurrentQ_create<Task_pack>(numParts);
//        timer = getCurrentTime();
//        for (int k = 0; k < LOOP_TIME; ++k) {
//            size_t offsetIdx = 0;
//            for (int i = 0; i < numParts; ++i) {
//                taskPack.workingLength = a->length / numParts;
//                taskPack.offset = offsetIdx;
//                offsetIdx += taskPack.workingLength;
//                if (a->length - offsetIdx < taskPack.workingLength || i == numParts - 1) {
//                    taskPack.workingLength = a->length - taskPack.offset;
//                    ConcurrentQ_enqueue(queue, taskPack);
//                    break;
//                }
//                ConcurrentQ_enqueue(queue, taskPack);
//            }
//            executionLayer(a, b, c, workCnt, s, e, queue);
//        }
//        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
//        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
//        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
//        printf("\n");
//
//        cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
//                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
//        clWaitForEvents(1, &clEvent);
//        ConcurrentQ_destroy(queue);
//    }
//
//    s = 0, e = 9;
//    cout << "Processors s: " << s << ", e: " << e << "\n";
//    printf("%-12s%-10s\t%s\n", "numParts", "Time", "Work count");
//    for (int j = 0; j < 14; ++j) {
//        memset(workCnt, 0, 20 * sizeof(float));
//        numParts = pow2(j);
//        ConcurrentQ<Task_pack> *queue = ConcurrentQ_create<Task_pack>(numParts);
//        timer = getCurrentTime();
//        for (int k = 0; k < LOOP_TIME; ++k) {
//            size_t offsetIdx = 0;
//            for (int i = 0; i < numParts; ++i) {
//                taskPack.workingLength = a->length / numParts;
//                taskPack.offset = offsetIdx;
//                offsetIdx += taskPack.workingLength;
//                if (a->length - offsetIdx < taskPack.workingLength || i == numParts - 1) {
//                    taskPack.workingLength = a->length - taskPack.offset;
//                    ConcurrentQ_enqueue(queue, taskPack);
//                    break;
//                }
//                ConcurrentQ_enqueue(queue, taskPack);
//            }
//            executionLayer(a, b, c, workCnt, s, e, queue);
//        }
//        printf("%-12d%-10f\t", numParts, (getCurrentTime() - timer) / LOOP_TIME);
//        for (int l = 0; l < 12; ++l) printf("%-5.2f ", workCnt[l] * 100 / (LOOP_TIME * numParts));
//        printf("\t diff: %f", calVecDiff(c->hostPtr, c_cpu, c->length));
//        printf("\n");
//
//        cl_wrapper::checkError(clEnqueueFillBuffer(exeCmdQueue, c->deviceBuffer, &cl_zero, sizeof(cl_int),
//                                                   0, c->length * sizeof(cl_int), 0, nullptr, &clEvent));
//        clWaitForEvents(1, &clEvent);
//        ConcurrentQ_destroy(queue);
//    }


    free(workCnt);
    cout << "Done! \n";
}


int main(int argc, char **argv) {
    const size_t w = 1024;
    const size_t h = 1024;
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

    cout << "Total time: " << timer << "s\n";//, diff: " << calVecDiff(c.hostPtr, c_cpu, dataLength) << "\n";


    return 0;
}

//    check scheduling queue, put in function shceduleLayer
//    std::cout << "====Data length: " << a->length << ", task queue====\n";
//    for (int i = 0; i < queue->size; ++i) {
//        std::cout << i << "\t" << queue->array[(queue->front + i) % queue->capacity].offset << "\t";
//        std::cout << queue->array[(queue->front + i) % queue->capacity].workingLength << "\t";
//        std::cout << queue->array[(queue->front + i) % queue->capacity].workingLength +
//                     queue->array[(queue->front + i) % queue->capacity].offset << "\n";
//    }
//    std::cout << "====(size " << queue->size;
//    std::cout << ", front " << queue->front;
//    std::cout << ", rear " << queue->rear;
//    std::cout << ")====\n\n";

// Code for testing Concurrent queue
//ConcurrentQ<float> *queue = ConcurrentQ_create<float>(5);
//ConcurrentQ_enqueue(queue, 0.0f);
//ConcurrentQ_enqueue(queue, 1.0f);
//ConcurrentQ_enqueue(queue, 2.0f);
//ConcurrentQ_print(queue);
//ConcurrentQ_enqueue(queue, 3.0f);
//ConcurrentQ_enqueue(queue, 4.0f);
//ConcurrentQ_enqueue(queue, 5.0f);
//ConcurrentQ_print(queue);
//float item;
//
//while (!ConcurrentQ_isEmpty(queue)) {
//ConcurrentQ_dequeueF(queue, &item);
//std::cout << "dequeueF: " << item << std::endl;
//ConcurrentQ_print(queue);
//}
//ConcurrentQ_enqueue(queue, 6.0f);
//ConcurrentQ_enqueue(queue, 7.0f);
//ConcurrentQ_enqueue(queue, 8.0f);
//ConcurrentQ_enqueue(queue, 9.0f);
//ConcurrentQ_enqueue(queue, 10.0f);
//ConcurrentQ_print(queue);
//while (!ConcurrentQ_isEmpty(queue)) {
//ConcurrentQ_dequeueR(queue, &item);
//std::cout << "dequeueR: " << item << std::endl;
//ConcurrentQ_print(queue);
//}
//
//for (int i = 0; i < queue->capacity; ++i) {
//std::cout << queue->array[i] << " ";
//}
//std::cout << "\n";
//
//ConcurrentQ_destroy(queue);
