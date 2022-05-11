//
// Created by pfxu on 3/13/19.
//
#include <iostream>
#include <cstdio>
#include <string.h>
#include <pthread.h>
#include <stdarg.h>
#include <sched.h>
#include <cassert>
#include <vector>

#include "CPUExecutor.h"
#include "utils.h"

int getNumCPUs() {
    long nprocs = -1;
    nprocs = sysconf(_SC_NPROCESSORS_ONLN);
    if (nprocs < 1) {
        flushed_printf("Could not determine the number of CPUs on line\n");
        flushed_printf("Error %d: %s\n", errno, strerror(errno));
    }

    return static_cast<int>(nprocs);
}

int setCurThreadAffinity(const std::vector<int>& cpus) {
    nice(-20);
    int policy = 0;
    struct sched_param param{};
    pthread_getschedparam(pthread_self(), &policy, &param);
    param.sched_priority = sched_get_priority_max(policy);
    pthread_setschedparam(pthread_self(), policy, &param);

    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);

    for (int cpu_id : cpus) {
        int numCPUs = getNumCPUs();
        if (cpu_id >= numCPUs) {
            flushed_printf("Fail to set CPU affinity...\n");
            flushed_printf("cpu_id (%d) doesn't exits, number of CPUs: %d\n",
                           cpu_id, numCPUs);
        }
        CPU_SET(cpu_id, &cpu_set);
    }
    pid_t thread_id = gettid();
    int err;
    if ((err = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set)) != 0) {
        //flushed_printf("Setting thread %d affinity on CPU %d failed %d\n", thread_id, cpu_id, err);
        //flushed_printf("Error %d: %s\n", errno, strerror(errno));
        return -1;
    }
    //long err = syscall(__NR_sched_setaffinity, thread_id, sizeof(cpu_set_t), &cpu_set);
    //if (err) printf("Setting affinity returns error(%ld): %s\n", err, strerror(errno));
    return 0;
}

int setCurThreadAffinity(int cpu_id) {
    // higher the priority
    nice(-20);
    int policy = 0;
    struct sched_param param{};
    pthread_getschedparam(pthread_self(), &policy, &param);
    param.sched_priority = sched_get_priority_max(policy);
    pthread_setschedparam(pthread_self(), policy, &param);

    int numCPUs = getNumCPUs();
    if (cpu_id >= numCPUs) {
        flushed_printf("Fail to set CPU affinity...\n");
        flushed_printf("cpu_id (%d) doesn't exits, number of CPUs: %d\n",
                       cpu_id, numCPUs);
    }

    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu_id, &cpu_set);

    pid_t thread_id = gettid();
    int err;
    if ((err = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set)) != 0) {
        //flushed_printf("Setting thread %d affinity on CPU %d failed %d\n", thread_id, cpu_id, err);
        //flushed_printf("Error %d: %s\n", errno, strerror(errno));
        return -1;
    }
    //long err = syscall(__NR_sched_setaffinity, thread_id, sizeof(cpu_set_t), &cpu_set);
    //if (err) printf("Setting affinity returns error(%ld): %s\n", err, strerror(errno));

    return 0;
}

cpu_set_t getCurThreadAffinity(pid_t threadID) {
    cpu_set_t cpuSet;

    if (sched_getaffinity(threadID, sizeof(cpu_set_t), &cpuSet) == -1) {
        perror("sched_getaffinity failed!\n");
        assert(false);
    } else {
        return cpuSet;
    }
}

void checkCurThreadAffinity(int thread_id) {
    cpu_set_t cpuSet = getCurThreadAffinity(gettid());
    int nproc, i;
    nproc = getNumCPUs();
    flushed_printf("Total number of CPUs = %d, Thread %d sched_getaffinity = ", nproc, thread_id);
    for (i = 0; i < nproc; i++) {
        if (CPU_ISSET(i, &cpuSet)) flushed_printf(" %d ", i);
    }
    flushed_printf("\n");
}

void *executeOnPthread(void *ptr) {
    cpu_task_pack *cpuTaskPack = (cpu_task_pack *) ptr;
    flushed_printf("Thread %d started.. \n", cpuTaskPack->thread_id);

    while (setCurThreadAffinity(cpuTaskPack->cpu_id) != 0) {
        //flushed_printf("Trying set thread %d affinity to %d ...\n",
        //               cpuTaskPack->thread_id, cpuTaskPack->cpu_id);
    }

    if (setCurThreadAffinity(cpuTaskPack->cpu_id)) {
        pthread_exit(nullptr);
    }
    checkCurThreadAffinity(cpuTaskPack->thread_id);

    cpuTaskPack->cpuKernel(cpuTaskPack->input,
                           cpuTaskPack->output,
                           cpuTaskPack->taskPack->offset,
                           cpuTaskPack->taskPack->length);

    flushed_printf("Thread %d Job Done!\n\n", cpuTaskPack->thread_id);
    pthread_exit(nullptr);
}

void executeOnCPU(pthread_t *worker, cpu_task_pack *cpuTaskPack) {
    if (pthread_create(worker, nullptr, &executeOnPthread, cpuTaskPack)) {
        flushed_printf("Pthread %d create failed...\n", cpuTaskPack->thread_id);
        flushed_printf("Error %d: %s\n", errno, strerror(errno));
    }
}
