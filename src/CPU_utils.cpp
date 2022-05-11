#include <tar.h>
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include "CPU_utils.h"
#include "utils.h"
#include <unistd.h>


int getNumCPUs() {
	long nprocs;
	nprocs = sysconf(_SC_NPROCESSORS_ONLN);
	if (nprocs < 1) {
		flushed_printf("Could not determine the number of CPUs on line\n");
		flushed_printf("Error %d: %s\n", errno, strerror(errno));
	}
	return static_cast<int>(nprocs);
}

int setThreadAffinity(u_int16_t cpuSet) {
	nice(-20);
	int policy = 0;
	struct sched_param param{};
	pthread_getschedparam(pthread_self(), &policy, &param);
	param.sched_priority = sched_get_priority_max(policy);
	pthread_setschedparam(pthread_self(), policy, &param);

	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);

	for (int i = 0; i < getNumCPUs(); i++) {
		if (IS_BIT_SET(cpuSet, i)) {
			CPU_SET(i, &cpu_set);
		}
	}
	pid_t thread_id = gettid();
	__unused int err;
	if ((err = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set)) != 0) {
//	EnableAllCores_();
//	flushed_printf("Setting thread %d affinity on CPU %d failed %d\n", thread_id, cpuSet, err);
//	flushed_printf("Error %d: %s\n", errno, strerror(errno));
		return -1;
	}
	//long err = syscall(__NR_sched_setaffinity, thread_id, sizeof(cpu_set_t), &cpu_set);
	//if (err) printf("Setting affinity returns error(%ld): %s\n", err, strerror(errno));
	return 0;
}

int setThreadAffinity(const std::vector<int> &cpus) {
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
			flushed_printf("cpu_id (%d) doesn't exits, number of CPUs: %d\n", cpu_id, numCPUs);
		}
		CPU_SET(cpu_id, &cpu_set);
	}
	pid_t thread_id = gettid();
	__unused int err;
	if ((err = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set)) != 0) {
//	EnableAllCores_();
		//flushed_printf("Setting thread %d affinity on CPU %d failed %d\n", thread_id, cpu_id, err);
		//flushed_printf("Error %d: %s\n", errno, strerror(errno));
		return -1;
	}
	//long err = syscall(__NR_sched_setaffinity, thread_id, sizeof(cpu_set_t), &cpu_set);
	//if (err) printf("Setting affinity returns error(%ld): %s\n", err, strerror(errno));
	return 0;
}

int setThreadAffinity(int cpu_id) {
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
	__unused int err;
	if ((err = sched_setaffinity(thread_id, sizeof(cpu_set), &cpu_set)) != 0) {
//	EnableAllCores_();
		//flushed_printf("Setting thread %d affinity on CPU %d failed %d\n", thread_id, cpu_id, err);
		//flushed_printf("Error %d: %s\n", errno, strerror(errno));
		return -1;
	}
	//long err = syscall(__NR_sched_setaffinity, thread_id, sizeof(cpu_set_t), &cpu_set);
	//if (err) printf("Setting affinity returns error(%ld): %s\n", err, strerror(errno));

	return 0;
}

cpu_set_t getThreadAffinity(pid_t threadID) {
	cpu_set_t cpuSet;
	if (sched_getaffinity(threadID, sizeof(cpu_set_t), &cpuSet) == -1) {
		perror("sched_getaffinity failed!\n");
		assert(false);
	} else {
		return cpuSet;
	}
}

__unused void printThreadAffinity(const char *dec) {
	cpu_set_t cpuSet = getThreadAffinity(gettid());
	int nproc, i;
	nproc = getNumCPUs();
	flushed_printf("%s thread affinity= ", dec);
	for (i = 0; i < nproc; i++) {
		if (CPU_ISSET(i, &cpuSet)) flushed_printf("%d ", i);
	}
	flushed_printf("\n");
}

void printThreadAffinity(int thread_id) {
	cpu_set_t cpuSet = getThreadAffinity(gettid());
	int nproc, i;
	nproc = getNumCPUs();
	flushed_printf("Thread %d sched_getaffinity = ", thread_id);
	for (i = 0; i < nproc; i++) {
		if (CPU_ISSET(i, &cpuSet)) flushed_printf(" %d ", i);
	}
	flushed_printf("\n");
}

