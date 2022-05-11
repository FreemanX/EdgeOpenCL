#pragma ide diagnostic ignored "readability-static-accessed-through-instance"
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"

#include <EDCL_Scheduler.h>
#include "Stg_Sequential.h"
#include "polybench.h"
#include "Stg_CA2020.h"

// For old version compatibility
#define globalSize getGlobalSize()
#define localSize getLocalSize()
#define offset getOffset()

//#define POLY syr
//#define INIT_POLY new POLY(edcl)
//#define timeL4 32.0
//#define timeB4 5.1
//#define timeG1 12.6

#define POLY GEMM
#define INIT_POLY new POLY(edcl)
#define timeL4 0.731703
#define timeB4 0.337742
#define timeG1 0.035099

//#define POLY mvt
//#define INIT_POLY new POLY(edcl)
//#define timeL4 1.386604
//#define timeB4 0.151086
//#define timeG1 0.091902

class CA2020_UnitTest {

 public:
	EDCL *edcl;
	CA2020 *ca2020;
	POLY *poly;
	ResourceManager rm;
	AtomicKernelSet polyAKS;

	CA2020_UnitTest() {
		// init
		edcl = new EDCL();
		rm = new ResourceManager_();
		ca2020 = new CA2020(edcl, rm);
		poly = INIT_POLY;
		polyAKS = poly->createBenchAKS();
	}

	~CA2020_UnitTest() { delete edcl; }

	void recalCK() const {
		BenchKernel *bench_kernels[13];
		bench_kernels[0] = new Conv3D(edcl);
		bench_kernels[1] = new Conv2D(edcl);
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
		std::unordered_map<std::string, EDKernel> kernelMap;
		std::vector<AtomicKernelSet> AKSs;
		for (auto &b : bench_kernels) { AKSs.push_back(b->createBenchAKS()); }
		for (auto &aks: AKSs) {
			for (auto &k : aks->kernels) {
				if (kernelMap.count(ca2020->genKernelID(k)) > 0) continue;
				kernelMap[ca2020->genKernelID(k)] = k;
			}
		}
		TDP_MAP polyTDP;
		ca2020->readTDPFromCSVFile("PolyTDP_New_distill.csv", polyTDP);
		for (auto &p : polyTDP) {
			for (auto &dp : p.second) {
				int totalCK = 0;
				for (auto &ck : dp.CK) {
					totalCK += ck;
				}
				if (totalCK > 0) continue;
				EDKernel kernel = kernelMap[dp.DP_ID];
				AtomicKernelSet chunks = edcl->slicingKernelPow2_extra_equal(kernel, dp.NC);
				int largestDim = 0;
				int lastUsedDevice = 0;
				size_t totalNumWG;
				{
					auto k = chunks->kernels[0];
					for (int kJ = 0; kJ < k->workDim; ++kJ) {
						if (k->globalSize[kJ] / k->localSize[kJ] > k->globalSize[largestDim] / k->localSize[largestDim])
							largestDim = kJ;
					}
					totalNumWG = k->globalSize[largestDim] / k->localSize[largestDim];
					for (int kI = 0; kI < 3; ++kI) {
						if (dp.CK_ratio[kI] > 0) lastUsedDevice = kI;
					}
				}
				std::unique_ptr<size_t[]> k_offset = std::make_unique<size_t[]>(kernel->workDim);
				std::unique_ptr<size_t[]> k_global = std::make_unique<size_t[]>(kernel->workDim);
				for (auto k : chunks->kernels) {
					memset(k_offset.get(), 0, k->workDim * sizeof(size_t));
					memset(k_global.get(), 0, k->workDim * sizeof(size_t));
					for (int kI = 0; kI < 3; ++kI) {
						if (dp.CK_ratio[kI] > 0) {
							EDKernel processorSlice = nullptr;
							for (int kJ = 0; kJ < k->workDim; ++kJ) { //kJ: loop through workDim
								if (kJ == largestDim) {
									size_t numSliceWG = floor((float)totalNumWG * dp.CK_ratio[kI]);
									if (numSliceWG > 0 || (kI == lastUsedDevice && k_offset[kJ] < k->globalSize[kJ])) {
										processorSlice = edcl->createKernelCopy(k);
										k_global[kJ] = numSliceWG * k->localSize[kJ];
									}
								} else k_global[kJ] = k->globalSize[kJ];
							}
							if (processorSlice == nullptr) continue;
							for (int kJ = 0; kJ < k->workDim; ++kJ) { // loop through work-dim, check if all WIs are included
								processorSlice->offset[kJ] = k_offset[kJ] + k->getOffset()[kJ];
								if (kI == lastUsedDevice && k->globalSize[kJ] - k_offset[kJ] > 0)
									processorSlice->globalSize[kJ] = k->globalSize[kJ] - k_offset[kJ];
								else
									processorSlice->globalSize[kJ] = k_global[kJ];
								if (kJ == largestDim) {
									k_offset[kJ] += k_global[kJ];
								}
							}
							dp.CK[kI] = processorSlice->globalSize[largestDim] / processorSlice->localSize[largestDim];
							assert(dp.CK[kI] > 0);
							delete processorSlice;
						}
					}
					delete k;
				}
				ca2020->printDP(dp);
			}
		}
		ca2020->writeTDPToCSVFile("PolyTDP_New_distill_ck.csv", polyTDP);
	}

//	void recalEC() const {
//		TDP_MAP polyTDP;
//		ca2020->readTDPFromCSVFile(TDP_FILE_NAME, polyTDP);
//		float i0 = 0.1167213;
//		for (auto &p : polyTDP) {
//			for (auto &dp: p.second) {
//				float t = 1 / dp.Prf;
//				float fl = CHOOSE(dp.nl >= 0, (float)(rm->getHeteroCPUFrequencies()[0].at(dp.fl)) / 1000000.0f, 0);
//				float fb = CHOOSE(dp.nb >= 0, (float)(rm->getHeteroCPUFrequencies()[1].at(dp.fb)) / 1000000.0f, 0);
//				float fg = CHOOSE(dp.ng >= 0, (float)(rm->getGPUFrequencies().at(dp.fg)) / 1000.0f, 0);
//				float tmp = (fl + fb * 2 + fg * 3) * t;
//				float I = sqrt((dp.EC / tmp) + i0 * i0);
//				float oldEC = dp.EC;
//				ca2020->recalEC(I, t, dp);
//				debug_printf(__func__, "Old EC: %f, I: %f, new EC: %f\n", oldEC, I, dp.EC);
//			}
//		}
//		ca2020->distillTDP("PolyTDP_Distilled_test.csv", polyTDP);
//	}

	void testBenchmarkKernelSlices() const {
		debug_printf(__func__, "Start...\n");
		for (auto &k : polyAKS->kernels) {
			DP dp{};
			dp.DP_ID = CA2020::genKernelID(k);
			dp.nl = -2;
			dp.nb = 2;
			dp.ng = 0;
			std::vector<EDQueue> queues{nullptr, nullptr, nullptr};
			if (dp.nl >= 0) queues[0] = (rm->createQueueLittleCore(dp.nl, edcl));
			if (dp.nb >= 0) queues[1] = (rm->createQueueBigCore(dp.nb, edcl));
			if (dp.ng >= 0) queues[2] = (edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue));
			std::vector<float> prfs{0, 0, 0};
			if (dp.nl >= 0) prfs[0] = 1.0 / timeL4;
			if (dp.nb >= 0) prfs[1] = 1.0 / timeB4;
			if (dp.ng >= 0) prfs[2] = 1.0 / timeG1;
			float sumPrfs = 0;
			for (int kI = 0; kI < 3; ++kI) { sumPrfs += prfs[kI]; }
			for (int i = 0; i < 3; ++i) { dp.CK_ratio[i] = prfs[i] * (1.0 / sumPrfs); }
			ca2020->findNC(k, dp, queues);
			auto aks = edcl->slicingKernelPow2_extra_equal(k, dp.NC);
			ca2020->setDP(dp, aks, queues);
			for (auto &q: queues) {
				if (q != nullptr)
					edcl->releaseEDQueue(q);
			}
			CA2020::printDPVertically(dp, rm);
		}
	}

	void testSetDP() const {
		debug_printf(__func__, "Start...\n");
		for (auto &k : polyAKS->kernels) {
			DP dp{};
			dp.DP_ID = CA2020::genKernelID(k);
			dp.nl = -1;
			dp.nb = -1;
			dp.ng = 0;
			std::vector<EDQueue> queues{nullptr, nullptr, nullptr};
			queues[0] = CHOOSE(dp.nl >= 0, rm->createQueueLittleCore(dp.nl, edcl), nullptr);
			queues[1] = CHOOSE(dp.nb >= 0, rm->createQueueBigCore(dp.nb, edcl), nullptr);
			queues[2] = CHOOSE(dp.ng >= 0, edcl->createDeviceCmdQueueProfilingEnabled(GPUQueue), nullptr);
			std::vector<float> prfs{0, 0, 0};
			if (dp.nl >= 0) prfs[0] = 1.0 / timeL4;
			if (dp.nb >= 0) prfs[1] = 1.0 / timeB4;
			if (dp.ng >= 0) prfs[2] = 1.0 / timeG1;

			float sumPrfs = 0;
			for (int kI = 0; kI < 3; ++kI) { sumPrfs += prfs[kI]; }
			for (int i = 0; i < 3; ++i) { dp.CK_ratio[i] = prfs[i] * (1.0 / sumPrfs); }

			auto freq_L = rm->getHeteroCPUFrequencies()[0];
			auto freq_B = rm->getHeteroCPUFrequencies()[1];
			auto freq_G = rm->getGPUFrequencies();
			int t_nfL = freq_L.size(); // total num little core frequency
			int t_nfB = freq_B.size();
			int t_nfG = freq_G.size();

			for (int nfG = 0; nfG < t_nfG; ++nfG) {
				for (int nfB = 0; nfB < t_nfB; ++nfB) {
					for (int nfL = 0; nfL < t_nfL; ++nfL) {
						dp.fl = nfL;
						dp.fb = nfB;
						dp.fg = nfG;
						rm->setCPUFrequency(0, dp.fl);
						rm->setCPUFrequency(rm->getNumCPU() - 1, dp.fb);
						rm->setGPUFrequency(dp.fg);
						ca2020->setDP(dp, polyAKS, queues);
						CA2020::printDPHead();
						CA2020::printDP(dp);
//			ca2020->printDPVertically(dp, rm);
					}
				}
			}
		}
	}

	void test_genE0() const {
		debug_printf(__func__, "Start...\n");
		ca2020->genI0();
		CA2020::printDPHead();
		for (auto &dp: ca2020->I0_DPV) {
			CA2020::printDP(dp);
		}
	}

	void test_read_write_DPCSVFile() const {
		debug_printf(__func__, "Start...\n");
		std::string testFile = "I0.csv";

		std::unordered_map<std::string, std::vector<DP>> TDP;
		if (ca2020->readTDPFromCSVFile(testFile.c_str(), TDP) == 0) {
//	  CA2020::printDPHead();
			for (auto &dp_p : TDP) {
				auto dp_v = dp_p.second;
				for (auto &dp : dp_v) {
					CA2020::printDP(dp);
				}
			}
		} else {
			err_printf(__func__, __LINE__, "CSV file doesn't exist!\n");
		}
		ca2020->writeTDPToCSVFile("I0_out.csv", TDP);
	}

	void test_genDPs() const {
		size_t SHRINK_FACTOR = 4;
		debug_printf(__func__, "Start...\n");
		std::unordered_map<std::string, std::vector<DP>> TDP;
		std::string TDPFile = "TDP_" + std::string(poly->getName()) + ".csv";
		ca2020->readTDPFromCSVFile(TDPFile.c_str(), TDP);
		for (auto &k : polyAKS->kernels) {
			if (SHRINK_FACTOR > 0) {
				for (int i = 0; i < k->workDim; ++i) {
					k->globalSize[i] /= SHRINK_FACTOR;
					k->localSize[i] /= 4;
				}
			}
			std::string kernelID = ca2020->genKernelID(k);
			if (TDP.count(kernelID) > 0) {
				ca2020->genDPs(k, TDP[kernelID], TDPFile.c_str());
			} else {
				std::vector<DP> DPs;
				ca2020->genDPs(k, DPs, TDPFile.c_str());
				TDP[kernelID] = DPs;
			}
		}
//	ca2020->printTDP(TDP);

//	time_t t = time(nullptr);
//	struct tm tm = *localtime(&t);
//	std::string filename = "PolyTDP_"
//		+ std::to_string(tm.tm_hour)
//		+ "_"
//		+ std::to_string(tm.tm_min)
//		+ "_"
//		+ std::to_string(tm.tm_sec)
//		+ ".csv";
//	ca2020->writeTDPToCSVFile(filename.c_str(), TDP);
	}

	void genAllDPs() const {
		BenchKernel *bench_kernels[13];
		bench_kernels[0] = new Conv3D(edcl);
		bench_kernels[1] = new Conv2D(edcl);
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
		TDP_MAP TDP;
		std::string TDPFile = "PolyTDP_New_1022_1200.csv";
		if (ca2020->readTDPFromCSVFile(TDPFile.c_str(), TDP) != 0)
			ca2020->writeTDPToCSVFile(TDPFile.c_str(), TDP);
		for (auto &b : bench_kernels) {
			debug_printf(__func__, "Generating DPs for %s\n", b->getName());
			AtomicKernelSet aks = b->createBenchAKS();
			for (auto &k : aks->kernels) {
				std::string kernelID = ca2020->genKernelID(k);
				if (TDP.count(kernelID) > 0) {
					ca2020->genDPs(k, TDP[kernelID], TDPFile.c_str());
				} else {
					std::vector<DP> DPs;
					ca2020->genDPs(k, DPs, TDPFile.c_str());
					TDP[kernelID] = DPs;
				}
			}
		}
		ca2020->distillTDP("PolyTDP_New_Distilled.csv", TDP);
	}

	void test_distillTDP() const {
		TDP_MAP gemmTDP;
		ca2020->readTDPFromCSVFile("TDP_gemm.csv", gemmTDP);
		ca2020->printTDP(gemmTDP);
		int wait;
		std::cin >> wait;
		ca2020->distillTDP("TDP_gemm_distilled.csv", gemmTDP);
	}

	void test_plan() const {
		EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
		Conv2D conv2d(edcl);
		Conv3D conv3d(edcl);
		MM2 mm2(edcl);
		ATAX atax(edcl);
		BICG bicg(edcl);
		Correlation corr(edcl);
		Covariance covar(edcl);
		FDTD2D fdtd2d(edcl);
		GEMM gemm(edcl);
		Gesummv gesummv(edcl);
		GramSchmidt gram_schmidt(edcl);
		mvt mvt(edcl);
		syr syrk(edcl);
		scheduler->submitAtomicKernelSets(
				conv2d.createBenchAKS(), conv3d.createBenchAKS(), mm2.createBenchAKS(),
				atax.createBenchAKS(), bicg.createBenchAKS(), corr.createBenchAKS(),
				covar.createBenchAKS(), fdtd2d.createBenchAKS(), gemm.createBenchAKS(),
				gesummv.createBenchAKS(), gram_schmidt.createBenchAKS(),
				mvt.createBenchAKS(), syrk.createBenchAKS());
		ca2020->plan(scheduler->atomicKernelSets, rm);
		debug_printf(__func__, "Planing time: %lf\n", ca2020->planingTime);
	}

	void test_genAllStageCDP() const {
		EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
		Conv2D conv2d(edcl);
		Gesummv gesummv(edcl);
		mvt mvt(edcl);
		syr syr(edcl);
		scheduler->submitAtomicKernelSets(conv2d.createBenchAKS(),
																			gesummv.createBenchAKS(),
																			mvt.createBenchAKS(),
																			syr.createBenchAKS());
		debug_printf(__func__, "Num Kernels: %d\n", scheduler->getNumEDKernels());

		std::vector<AtomicKernelSet> priorityStages;
		scheduler->InitSchedulerEnv(priorityStages);
		ca2020->plan(priorityStages, rm);

		int stageCnd = 0;
		for (auto stage: priorityStages) {
			debug_printf(__func__, "Stage %d\n", stageCnd++);
			std::vector<std::vector<DP> > stage_CDPs;
			double time = getCurrentTime();
			ca2020->genAllStageCDP(stage_CDPs, stage, 8);
			time = getCurrentTime() - time;
			debug_printf("genAllStageCDP", "Num CDPs: %d, time: %lf\n", stage_CDPs.size(), time);
			auto start = 3 * (int)stage_CDPs.size() / 4;
			for (int i = start; i < stage_CDPs.size() && i < 5 + start; i++) {
				debug_printf(__func__, "CDP[%d]...\n", i);
				ca2020->printCDP(stage_CDPs[i]);
			}
			std::sort(stage_CDPs.begin(), stage_CDPs.end(), PrfIsBetter);
			for (int i = 0; i < 5; i++) {
				debug_printf(__func__, "Prf sorted CDP[%d]...\n", i);
				ca2020->printCDP(stage_CDPs[i]);
			}
			std::sort(stage_CDPs.begin(), stage_CDPs.end(), ECIsBetter);
			for (int i = 0; i < 5; i++) {
				debug_printf(__func__, "EC sorted CDP[%d]...\n", i);
				ca2020->printCDP(stage_CDPs[i]);
			}
			std::sort(stage_CDPs.begin(), stage_CDPs.end(), ProductIsBetter);
			for (int i = 0; i < 5; i++) {
				debug_printf(__func__, "Product sorted CDP[%d]...\n", i);
				ca2020->printCDP(stage_CDPs[i]);
			}
		}
	}

	void compare_getCDP() const {
		EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
		Conv2D conv2d(edcl);
		Gesummv gesummv(edcl);
		mvt mvt(edcl);
		syr syr(edcl);
		scheduler->submitAtomicKernelSets(conv2d.createBenchAKS(),
																			gesummv.createBenchAKS(),
																			mvt.createBenchAKS(),
																			syr.createBenchAKS());
		debug_printf(__func__, "Num Kernels: %d\n", scheduler->getNumEDKernels());

		std::vector<AtomicKernelSet> priorityStages;
		scheduler->InitSchedulerEnv(priorityStages);
		ca2020->plan(priorityStages, rm);

		int stageCnd = 0;
		for (auto stage: priorityStages) {
			debug_printf(__func__, "Stage %d\n", stageCnd++);

			double time = getCurrentTime();
			std::vector<std::vector<DP> > stage_CDPs;
			ca2020->genAllStageCDP(stage_CDPs, stage, 8);
			std::sort(stage_CDPs.begin(), stage_CDPs.end(), PrfIsBetter);
			time = getCurrentTime() - time;
			debug_printf("genAllStageCDP", "Num CDPs: %d, time: %lf\n", stage_CDPs.size(), time);
			for (int i = 0; i < 5; i++) {
				debug_printf(__func__, "Prf sorted CDP[%d]...\n", i);
				ca2020->printCDP(stage_CDPs[i]);
			}

//	  time = getCurrentTime();
//	  std::vector<std::vector<DP> > CDPTopN;
//	  int n = 5;
//	  ca2020->genBestStageCDP_V2(n, CDPTopN, stage, sortByPrf);
//	  time = getCurrentTime() - time;
//	  debug_printf("genBestStageCDP_v2", "Num CDPs: %d, time: %lf\n", CDPTopN.size(), time);
//	  for (int i = 0; i < CDPTopN.size(); i++) {
//		debug_printf(__func__, "Prf top %d CDP...\n", i);
//		printCDP(CDPTopN[i]);
//	  }

			time = getCurrentTime();
			std::vector<DP> bestCDP;
			ca2020->genBestStageCDP(bestCDP, stage, 8);
			time = getCurrentTime() - time;
			debug_printf("genBestStageCDP", "Time: %lf\n", time);
			debug_printf(__func__, "Prf sorted CDP...\n");
			ca2020->printCDP(bestCDP);

		}
	}

	void test_execute() const {
		EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
		Conv2D conv2d(edcl);
		Gesummv gesummv(edcl);
		mvt mvt(edcl);
		syr syr(edcl);
		GramSchmidt gram_schmidt(edcl);
		FDTD2D fdtd2d(edcl);
		scheduler->submitAtomicKernelSets(
//				poly->createBenchAKS(),
//				syr.createBenchAKS(),
//				syr.createBenchAKS(),
//				gesummv.createBenchAKS(),
				conv2d.createBenchAKS()
//									  syr.createBenchAKS(),
//									  syr.createBenchAKS(),
//				gesummv.createBenchAKS(),
//				mvt.createBenchAKS(),
//				syr.createBenchAKS(),
//				gram_schmidt.createBenchAKS()
//				fdtd2d.createBenchAKS()
		);
		debug_printf(__func__, "Num Kernels: %d\n", scheduler->getNumEDKernels());

		SequentialImproved stg_seq(edcl);
		stg_seq.setExecutionDevice(GPU);
		scheduler->executeKernels(stg_seq);
		debug_printf(__func__, "Seq execution time: %lf\n", stg_seq.executionTime);

		std::vector<AtomicKernelSet> priorityStages;
		scheduler->InitSchedulerEnv(priorityStages);
		ca2020->setPreference(PERFORMANCE);
		ca2020->plan(priorityStages, rm);
		ca2020->execute(rm);
		debug_printf(__func__, "Seq execution time: %lf\n", stg_seq.executionTime);
		debug_printf(__func__, " CA2020 execution time: %lf\n", ca2020->executionTime);
	}

	void test_all() const {
		EDCL_Scheduler *scheduler = GET_SCHEDULER(edcl);
		Conv2D conv2d(edcl);
		Conv3D conv3d(edcl);
		MM2 mm2(edcl);
		ATAX atax(edcl);
		BICG bicg(edcl);
		Correlation corr(edcl);
		Covariance covar(edcl);
		FDTD2D fdtd2d(edcl);
		GEMM gemm(edcl);
		Gesummv gesummv(edcl);
		GramSchmidt gram_schmidt(edcl);
		mvt mvt(edcl);
		syr syrk(edcl);
		scheduler->submitAtomicKernelSets(
				conv3d.createBenchAKS(),
				conv2d.createBenchAKS(),
				mm2.createBenchAKS(),
				atax.createBenchAKS(), bicg.createBenchAKS(), corr.createBenchAKS(),
				covar.createBenchAKS(), fdtd2d.createBenchAKS(), gemm.createBenchAKS(),
				gesummv.createBenchAKS(), gram_schmidt.createBenchAKS(),
				mvt.createBenchAKS(), syrk.createBenchAKS());
		for (auto &aks : scheduler->atomicKernelSets) {
			for (auto &k : aks->kernels) {
				k->schedulePriority = 0;
			}
		}
		scheduler->executeKernels(*ca2020);
		debug_printf(__func__, " Planing time: %lf\n", ca2020->planingTime);
		debug_printf(__func__, " CA2020 execution time: %lf\n", ca2020->executionTime);
	}
};

int main() {
//  setThreadAffinity(0);
	CA2020_UnitTest unit_test;
	std::cout << "1: testSetDP()" << "\n";
	std::cout << "2: testBenchmarkKernelSlices()" << "\n";
	std::cout << "3: test_genE0()" << "\n";
	std::cout << "4: test_read_write_DPCSVFile()" << "\n";
	std::cout << "5: test_genDPs()" << "\n";
	std::cout << "6: genAllDPs()" << "\n";
	std::cout << "7: test_distillTDP()" << "\n";
	std::cout << "8: test_plan()" << "\n";
	std::cout << "9: genAllStageCDP()" << "\n";
	std::cout << "10: compare_getCDP()" << "\n";
	std::cout << "11: test_execute()" << "\n";
	std::cout << "12: test_all()" << "\n";

//  PRINT_RED_B
//  std::cout << "920: recalEC()" << "\n";
//  PRINT_DEFAULT

	int chose = 0;
	std::cin >> chose;
	switch (chose) {
		TEST_CASE(12, unit_test.test_all())
		TEST_CASE(11, unit_test.test_execute())
		TEST_CASE(10, unit_test.compare_getCDP())
		TEST_CASE(9, unit_test.test_genAllStageCDP())
//		TEST_CASE(920, unit_test.recalEC())
		TEST_CASE(930, unit_test.recalCK());
		case 1: {
			unit_test.testSetDP();
			break;
		}
		case 2: {
			unit_test.testBenchmarkKernelSlices();
			break;
		}
		case 3: {
			unit_test.test_genE0();
			break;
		}
		case 4: {
			unit_test.test_read_write_DPCSVFile();
			break;
		}
		case 5: {
			unit_test.test_genDPs();
			break;
		}
		case 6: {
			unit_test.genAllDPs();
			break;
		}
		case 7: {
			unit_test.test_distillTDP();
			break;
		}
		case 8: {
			unit_test.test_plan();
			break;
		}
		default: {
			err_printf("Unknown Option", chose, "\n");
			break;
		}

	}

	unit_test.rm->setSoCMaxFrequency();
	return 0;
}
