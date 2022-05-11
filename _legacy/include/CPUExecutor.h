//
// Created by pfxu on 3/13/19.
//

#ifndef PROJECT_CPUEXECUTOR_H
#define PROJECT_CPUEXECUTOR_H

#include "Executor.h"
#include "TaskPack.h"
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <sys/syscall.h>

//#ifndef CPU_ZERO
//#define CPU_SETSIZE 1024
//#define __NCPUBITS  (8 * sizeof (unsigned long))
//typedef struct
//{
//    unsigned long __bits[CPU_SETSIZE / __NCPUBITS];
//} cpu_set_t;
//
//#define CPU_SET(cpu, cpusetp) \
//  ((cpusetp)->__bits[(cpu)/__NCPUBITS] |= (1UL << ((cpu) % __NCPUBITS)))
//#define CPU_ZERO(cpusetp) \
//  memset((cpusetp), 0, sizeof(cpu_set_t))
//#elsesu -c netcfg rndis0 dhcp"
//#define CPU_SET(cpu, cpustep) ((void)0)
//#define CPU_ZERO(cpu, cpustep) ((void)0)
//#endif

int getNumCPUs();

int setCurThreadAffinity(int cpu_id);

int setCurThreadAffinity(const std::vector<int>& cpu_set);

void checkCurThreadAffinity(int thread_id);

void executeOnCPU(pthread_t *worker, cpu_task_pack *cpuTaskPack);

#endif //PROJECT_CPUEXECUTOR_H
