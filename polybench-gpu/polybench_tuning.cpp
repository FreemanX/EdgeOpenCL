#include "polybench.h"
#include <unordered_map>
#include <unistd.h>
#include "EDCL_ResourceManager.h"

BenchKernel *bench_kernels[13];
EDCL *edcl;

struct KernelBest {
	EDKernel kernel{};
	std::string device;
	size_t bestLocalSize[3] = {0, 0, 0};
	double minTime = MAXFLOAT;
	void printLocalSize() {
		std::cout << "[ ";
		for (int i = 0; i < kernel->workDim; ++i) {
			std::cout << bestLocalSize[i] << " ";
		}
		std::cout << "]";
	}
};

std::vector<KernelBest> KernelBestsB4;
std::vector<KernelBest> KernelBestsB2;
std::vector<KernelBest> KernelBestsB1;
std::vector<KernelBest> KernelBestsL4;
std::vector<KernelBest> KernelBestsGPU;

double exeKernel(EDQueue Queue, EDKernel k) {
	int loopNum = 10;
	double avgTime = 0;
	try {
		for (int i = 0; i < loopNum; ++i) {
			edcl->confirmExeEnvAndExecuteSingleKernel(Queue, k, 0, nullptr, &k->kernelEvent);
			k->waitForKernelEvent(Queue);
			avgTime += k->kernelTime;
		}
		avgTime = avgTime / loopNum;
	} catch (std::runtime_error &e) {
		avgTime = MAXFLOAT;
	}
	return avgTime;
}

// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()
void tuningLocalSize(EDQueue Queue, std::string &&Device, std::vector<KernelBest> &KernelBests) {
	double avgTime;
	for (auto &b : bench_kernels) {
		flushed_printf(" Running %s\n", b->getName());
		AtomicKernelSet aks = b->createBenchAKS();
		std::unordered_map<std::string, EDKernel> kernelRecord; // eliminate duplicated kernels
		for (auto &k : aks->kernels) {
			if (kernelRecord.count(k->name) > 0) continue;
			sleep(60);
			KernelBest kb;
			kb.kernel = k;
			kb.device = Device;

			if (k->workDim == 1) {
				for (int i = 0; i <= 10; ++i) {
					k->localSize[0] = pow2(i);

					avgTime = exeKernel(Queue, k);
					flushed_printf("\t%s single WorkDim(%d)[ ", k->name.c_str(), k->workDim);
					for (int kI = 0; kI < k->workDim; kI++) {
						flushed_printf("%d ", k->localSize[kI]);
					}
					if (avgTime == MAXFLOAT) {
						flushed_printf("]:\n");
					} else {
						flushed_printf("]:%lf\n", avgTime);
						if (avgTime < kb.minTime) {
							kb.minTime = avgTime;
							kb.bestLocalSize[0] = k->localSize[0];
						}
					}
				}
			} else if (k->workDim == 2) {
				for (int i = 1; i <= 10; ++i) {
					k->localSize[0] = pow2(i);
					for (int j = 1; j <= 10; ++j) {
						k->localSize[1] = pow2(j);

						avgTime = exeKernel(Queue, k);

						flushed_printf("\t%s single WorkDim(%d)[ ", k->name.c_str(), k->workDim);
						for (int kI = 0; kI < k->workDim; kI++) {
							flushed_printf("%d ", k->localSize[kI]);
						}
						if (avgTime == MAXFLOAT) {
							flushed_printf("]:\n");
						} else {
							flushed_printf("]:%lf\n", avgTime);
							if (avgTime < kb.minTime) {
								kb.minTime = avgTime;
								kb.bestLocalSize[0] = k->localSize[0];
								kb.bestLocalSize[1] = k->localSize[1];
							}
						}
					}
				}
			}
			kernelRecord[k->name] = k;
			KernelBests.push_back(kb);
		}
	}
}

void printBest(std::vector<KernelBest> &KernelBests) {
	std::cout << "===Kernel Best LWG vs. Time Report===\n";
	for (auto &kb : KernelBests) {
		EDKernel k = kb.kernel;
		std::cout << "\t" << k->name << " best run on " << kb.device << " WorkDim("
							<< k->workDim << ")[ ";
		for (int kI = 0; kI < k->workDim; kI++) {
			std::cout << kb.bestLocalSize[kI] << " ";
		}
		std::cout << "]" << ":" << kb.minTime << "\n";
	}
}

void printCSV() {
	// kernel, b4, b2, b1, l4, GPU
	std::cout << "===Kernel Best LWG CSV Report===\n";
	std::cout << "kernel,b4,b2,b1,l4,GPU\n";
	for (int i = 0; i < KernelBestsB4.size(); ++i) {
		auto k_name = KernelBestsB4[i].kernel->name;
		std::cout << k_name;
		std::cout << ",";
		KernelBestsB4[i].printLocalSize();
		std::cout << ",";
		KernelBestsB2[i].printLocalSize();
		std::cout << ",";
		KernelBestsB1[i].printLocalSize();
		std::cout << ",";
		KernelBestsL4[i].printLocalSize();
		std::cout << ",";
		KernelBestsGPU[i].printLocalSize();
		std::cout << "\n";
	}
}

int main() {
	edcl = new EDCL();
	bench_kernels[0] = new Conv2D(edcl);
	bench_kernels[1] = new Conv3D(edcl);
	bench_kernels[2] = new MM2(edcl);
	bench_kernels[3] = new ATAX(edcl);
	bench_kernels[4] = new BICG(edcl);
	bench_kernels[5] = new Correlation(edcl);
	bench_kernels[6] = new Covariance(edcl);
	bench_kernels[7] = new FDTD2D(edcl);
	bench_kernels[8] = new GEMM(edcl);
	bench_kernels[9] = new Gesummv(edcl);
	bench_kernels[10] = new GramSchmidt(edcl);
	bench_kernels[11] = new mvt(edcl);
	bench_kernels[12] = new syr(edcl);

	auto rm = new ResourceManager_();
	EDQueue Queue = rm->createUniqQueueBigCore(2, edcl);
	flushed_printf("Tuning CPU(Big 4)...\n");
	tuningLocalSize(Queue, "B4", KernelBestsB4);
	rm->releaseUniqQueue(Queue, edcl);

	Queue = rm->createUniqQueueBigCore(1, edcl);
	flushed_printf("Tuning CPU(Big 2)...\n");
	tuningLocalSize(Queue, "B2", KernelBestsB2);
	rm->releaseUniqQueue(Queue, edcl);

	Queue = rm->createUniqQueueBigCore(0, edcl);
	flushed_printf("Tuning CPU(Big 1)...\n");
	tuningLocalSize(Queue, "B1", KernelBestsB1);
	rm->releaseUniqQueue(Queue, edcl);

	Queue = rm->createUniqQueueLittleCore(2, edcl);
	flushed_printf("Tuning CPU(Little 4)...\n");
	tuningLocalSize(Queue, "L4", KernelBestsL4);
	rm->releaseUniqQueue(Queue, edcl);

	Queue = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
	flushed_printf("Tuning GPU...\n");
	tuningLocalSize(Queue, "GPU", KernelBestsGPU);

	printBest(KernelBestsB4);
	printBest(KernelBestsB2);
	printBest(KernelBestsB1);
	printBest(KernelBestsL4);
	printBest(KernelBestsGPU);

	printCSV();
	return 0;
}