#include <unistd.h>
#include "utils.h"
#include "CPU_utils.h"

int holdingTime = 300;

static void *holdCore(void *arg) {
	nice(39);
	int coreID = *(int *)arg;
	while (true) {
		int err = setThreadAffinity(coreID);
		if (err == 0) break;
	}
	debug_printf(__func__, "Holding core %d for %d sec...\n", coreID, holdingTime);
	nice(39);
	sleep(holdingTime);
	pthread_exit(nullptr);
}


// TODO lower priority

int main(int argc, char **argv) {
	nice(39);
	if (argc > 1) holdingTime = std::stoi(argv[1]);
	int bID[4] = {4, 5, 6, 7};
	pthread_t pthreads[40];
	for (int kI = 0; kI < 10; ++kI) {
		for (int j = 0; j < 4; ++j) {
			auto err = pthread_create(&pthreads[kI * 4 + j], nullptr, holdCore, &bID[j]);
			if (err) {
				err_printf("pthread_create", err, "pthread creation error!\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	for (long pthread : pthreads) {
		pthread_join(pthread, nullptr);
	}
	return 0;
}