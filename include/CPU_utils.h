#ifndef EDCL_SRC_CPU_UTILS_H_
#define EDCL_SRC_CPU_UTILS_H_
#include "utils.h"
#include <vector>

// TODO: get CPU available frequencies and power governor. Set CPU frequency


int getNumCPUs();

// returns 0 on success, -1 otherwise
int setThreadAffinity(u_int16_t cpuSet);

// returns 0 on success, -1 otherwise
// Redundant function...
int setThreadAffinity(int cpu_id);

// returns 0 on success, -1 otherwise
// Redundant function...
int setThreadAffinity(const std::vector<int> &cpus);

void printThreadAffinity(int thread_id);
__unused void printThreadAffinity(const char *dec);

#endif //EDCL_SRC_CPU_UTILS_H_
