#include "polybench.h"
#include "csv.h"
#include <unordered_map>
#include <unistd.h>
#include "EDCL_Strategy.h"
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
BenchKernel *bench_kernels[13];

EDCL *edcl;

double exeKernel(EDQueue Queue, EDKernel k, int loopNum) {
  double avgTime = 0;
  try {
	for (int i = 0; i < loopNum; ++i) {
	  edcl->confirmExeEnvAndExecuteSingleKernel(Queue, k, 0, nullptr, &k->kernelEvent);
	  k->waitForKernelEvent(Queue);
	  avgTime += k->kernelTime;
	}
	avgTime = avgTime / loopNum;
  } catch (std::runtime_error &e) {
	avgTime = 0;
  }
  return avgTime;
}

void print_map(std::unordered_map<std::string, std::vector<size_t>> const &m) {
  for (const auto &it : m) {
	std::cout << "{" << it.first << ": [";
	for (const auto &num : it.second)
	  std::cout << num << " ";
	std::cout << "]}\n";
  }
}

void addLocalWG(std::unordered_map<std::string, std::vector<size_t>> &map,
				std::string &KernelName,
				std::string &WGSizeStr) {
  if (map.count(KernelName) > 0) return;
  std::vector<size_t> WGSize;
  string2Size_t(WGSizeStr, WGSize);
  map[KernelName] = WGSize;
}

std::unordered_map<std::string, std::vector<size_t>> b4Map;
std::unordered_map<std::string, std::vector<size_t>> b2Map;
std::unordered_map<std::string, std::vector<size_t>> b1Map;
std::unordered_map<std::string, std::vector<size_t>> l4Map;
std::unordered_map<std::string, std::vector<size_t>> gpuMap;
int initKernelWGMap() {
  try {
	io::CSVReader<6> in("KernelWG.csv");
	in.read_header(io::ignore_extra_column, "kernel", "b4", "b2", "b1", "l4", "GPU");
	std::string kernel;
	std::string b4;
	std::string b2;
	std::string b1;
	std::string l4;
	std::string gpu;
	while (in.read_row(kernel, b4, b2, b1, l4, gpu)) {
	  addLocalWG(b4Map, kernel, b4);
	  addLocalWG(b2Map, kernel, b2);
	  addLocalWG(b1Map, kernel, b1);
	  addLocalWG(l4Map, kernel, l4);
	  addLocalWG(gpuMap, kernel, gpu);
	}
  } catch (io::error::can_not_open_file &e) {
	err_printf(__func__, __LINE__, "%s\n", e.what());
	return -1;
  }
  return 0;
}

void testKernelScalability(EDQueue Queue, std::unordered_map<std::string, std::vector<size_t>> &map) {
  double avgTime;
  for (auto &b : bench_kernels) {
	AtomicKernelSet aks = b->createBenchAKS();
	std::unordered_map<std::string, EDKernel> kernelRecord; // eliminate duplicated kernels
	for (auto &k : aks->kernels) {
	  if (kernelRecord.count(k->name) > 0) continue;
	  { // Dry run
		k->setLocalSize(map[k->name].data());
		edcl->confirmExeEnvAndExecuteSingleKernel(Queue, k, 0, nullptr, &k->kernelEvent);
		k->waitForKernelEvent(Queue);
	  }

	  int loopNum = MIN(MAX(10, floor(10.0 / k->kernelTime)), 200);
//	  sleep(20); // let the SoC cool down
	  kernelRecord[k->name] = k;
	  std::cout << " Running " << b->getName() << " global size" << k->name << "[ ";
	  for (int dim = 0; dim < k->workDim; ++dim) std::cout << k->getGlobalSize()[dim] << " ";
	  flushed_printf("], loopNum: %d\n", loopNum);

	  for (int kI = 0; kI < 11; ++kI) {
		AtomicKernelSet sub_kernelSet = edcl->slicingKernelPow2_extra_equal(k, kI);
		if (sub_kernelSet == nullptr) break;
		// time a chunk
		EDKernel testChunk = sub_kernelSet->kernels[0];
		testChunk->kernelTime = 0;
		avgTime = exeKernel(Queue, testChunk, loopNum);
		flushed_printf("\t");
		std::cout << testChunk->name << " single slice( ";
		for (int dim = 0; dim < k->workDim; ++dim) std::cout << testChunk->getGlobalSize()[dim] << " ";
		std::cout << ") time: " << avgTime << "\n";
		flushed_printf("");

		// time the kernel divided into pow2(kI) chunks
		avgTime = 0;
		for (int kJ = 0; kJ < loopNum; ++kJ) {
		  executeAtomicKernelSet(edcl, sub_kernelSet, Queue);
		  avgTime += sub_kernelSet->atomicExeTime;
		}
		avgTime = avgTime / loopNum;
		flushed_printf("\t");
		std::cout << k->name << " num slices: " << sub_kernelSet->getNumKernels() << " total time: " << avgTime << "\n";

//		delete sub_kernelSet;
		flushed_printf("");
	  }
	}
//	delete aks;
  }
}

int main() {
  int ret;
  if ((ret = initKernelWGMap()) != 0)
	return ret;
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
  EDQueue Queue;

  flushed_printf("Testing GPU...\n");
  Queue = edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue);
  testKernelScalability(Queue, gpuMap);

  Queue = edcl->createSubDeviceExeQueue(2, 1);
  flushed_printf("Testing CPU(Big 4)...\n");
  testKernelScalability(Queue, b4Map);

  Queue = edcl->createSubDeviceExeQueue(1, 3);
  flushed_printf("Testing CPU(Big 2)...\n");
  testKernelScalability(Queue, b2Map);

  Queue = edcl->createSubDeviceExeQueue(0, 7);
  flushed_printf("Testing CPU(Big 1)...\n");
  testKernelScalability(Queue, b1Map);

  Queue = edcl->createSubDeviceExeQueue(2, 0);
  flushed_printf("Testing CPU(Little 4)...\n");
  testKernelScalability(Queue, l4Map);

  return 0;
}