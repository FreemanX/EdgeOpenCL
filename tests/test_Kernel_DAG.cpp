#include "EDCL.h"
#include "EDCL_Scheduler.h"

const char *source =
#include "vecAdd.cl"
unsigned int vDim = 32;
uint n = vDim * vDim;
size_t nByte = n * sizeof(float);
size_t glob_VA[1] = {n};
size_t local_VA[1] = {1024};

#define FLOAT_PTR(EDBuffer) HOST_PTR(float, EDBuffer)
EDCL *edcl;
EDProgram prog;
EDBuffer bA;
EDBuffer bB;
EDBuffer bC;
EDCL_Scheduler *scheduler;
EDKernel genKernel(std::string &&name) {
  EDKernel vecAdd = edcl->createKernel(prog, "vecAdd", 4);
  vecAdd->configKernel(1, glob_VA, local_VA, bA, bB, bC, &n);
  vecAdd->name = name;
  return vecAdd;
}

int main() {
  edcl = new EDCL();
  scheduler = EDCL_Scheduler::getInstance(edcl);
  // Init Buffers
  bA = edcl->createBuffer(nByte);
  bB = edcl->createBuffer(nByte);
  bC = edcl->createBuffer(nByte);
  randArrayGenerator<float>(0.001, 1, (float *)bA->hostPtr, n);
  randArrayGenerator<float>(0.001, 1, (float *)bB->hostPtr, n);
  memset(bC->hostPtr, 0, nByte);
  prog = edcl->createProgram(&source);
  auto A = genKernel("A");
  auto B = genKernel("B");
  auto C = genKernel("C");
  auto D = genKernel("D");
  auto E = genKernel("E");
  auto F = genKernel("F");
  auto G = genKernel("G");
  auto H = genKernel("H");
  A->addSubscriber(edcl, B);
  A->addSubscriber(edcl, C);
  A->addSubscriber(edcl, D);
  B->addSubscriber(edcl, D);
  C->addSubscriber(edcl, D);
  C->addSubscriber(edcl, E);
  F->addSubscriber(edcl, E);
  F->addSubscriber(edcl, G);
  auto aks = edcl->createAtomicKernelSet(8);
  aks->addKernel(A, B, C, D, E, F, G, H);
  std::vector<int> path;
  for (auto &k : aks->kernels) {
	int suggested = k->schedulePriority - 1;
//	debug_printf("Visiting", "%s...\n", k->name.c_str());
	for (auto &subscriber : k->subscribers) {
//	  debug_printf("Visiting subscriber", "%s...\n", subscriber->name.c_str());
	  subscriber->updatePriority(suggested, path);
	  path.clear();
	}
  }

  // sort kernel based on priority in descending order
  std::sort(aks->kernels.begin(), aks->kernels.end(),
			[](const EDKernel &k1, const EDKernel &k2) -> bool {
			  return k1->schedulePriority > k2->schedulePriority;
			});
  for (auto &k : aks->kernels) {
	flushed_printf("%s: %d\n", k->name.c_str(), k->schedulePriority);
  }

  std::vector<AtomicKernelSet> priorityStage;
  int lastPriority = 1;
  for (auto &k : aks->kernels) {
	if (lastPriority != k->schedulePriority) {
	  priorityStage.push_back(edcl->createAtomicKernelSet(1));
	  lastPriority = k->schedulePriority;
	}
	priorityStage.back()->addKernel(k);
  }

  for (int i = 0; i < priorityStage.size(); i++) {
	flushed_printf("Stage %d:\n", i);
	for (auto &k : priorityStage[i]->kernels) {
	  flushed_printf("%s: %d\n", k->name.c_str(), k->schedulePriority);
	}
  }

  for (auto &k : aks->kernels) {
	delete k;
  }
  return 0;
}