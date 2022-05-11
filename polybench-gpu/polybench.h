#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#ifndef EDCL_POLYBENCH_GPU_POLYBENCH_H_
#define EDCL_POLYBENCH_GPU_POLYBENCH_H_
#include "EDCL.h"
#include "utils.h"

typedef float DATA_TYPE;
#define PTR(EDBuffer) HOST_PTR(DATA_TYPE, EDBuffer)
#define RT_NAME(name) const char *getName() override { return name; }
#define sqrt_of_array_cell(x, j) sqrt(x[j])

class BenchKernel {
 public:
  explicit BenchKernel(EDCL *edcl) {
	this->edcl_ = edcl;
  }
  virtual AtomicKernelSet createBenchAKS() = 0;
  virtual const char *getName() = 0;
  virtual double nativeKernel(int NumThreads) = 0;
  const char **getKernelSource() { return &this->kernelSource; }
 protected:
  EDCL *edcl_;
  EDProgram program_ = nullptr;
  std::string kernelName;
  uint numArgs{};
  uint workDim{};
  size_t globalWorkSize[3]{0, 0, 0};
  size_t localWorkSize[3]{0, 0, 0};
  const char *kernelSource{};
};

class syr : public BenchKernel {
 public:
  RT_NAME("syr")
  explicit syr(EDCL *edcl, int n = 2048, int m = 2048) : BenchKernel(edcl) {
	N = n;
	M = m;
	alpha = 1;
	beta = 1;
	A = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	B = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	C = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	init_arrays(PTR(A), PTR(B), PTR(C));
	numArgs = 7;
	workDim = 2;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
	kernelSource =
#include "Kernels/syr2k.cl"
	program_ = edcl->createProgram(&kernelSource);
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(1);
	EDKernel k = edcl_->createKernel(program_, "syr2k_kernel", 7);
	k->configKernel(workDim, globalWorkSize, localWorkSize,
					A, B, C, &alpha, &beta, &M, &N);
	aks->addKernel(k);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	DATA_TYPE *C_ = PTR(C);
	int i, j, k;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, C_, M, N, alpha, beta) private(i, j, k)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		  C_[i * N + j] *= beta;
		}
	  }
#pragma omp for schedule(static)
	  for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		  for (k = 0; k < M; k++) {
			C_[i * N + j] += alpha * A_[i * M + k] * B_[j * M + k];
			C_[i * N + j] += alpha * B_[i * M + k] * A_[j * M + k];
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  EDBuffer A, B, C;
  int M, N;
  DATA_TYPE alpha, beta;
  size_t DIM_LOCAL_WORK_GROUP_X = 32;
  size_t DIM_LOCAL_WORK_GROUP_Y = 8;
  void init_arrays(DATA_TYPE *A_, DATA_TYPE *B_, DATA_TYPE *C_) const {
	int i, j;
#pragma omp parallel for private(i, j) shared(A_, B_, C_, N) num_threads(4)
	for (i = 0; i < N; i++) {
	  for (j = 0; j < N; j++) {
		C_[i * N + j] = ((DATA_TYPE)i * (DATA_TYPE)j + 2) / (DATA_TYPE)N;
	  }
	  for (j = 0; j < M; j++) {
		A_[i * N + j] = ((DATA_TYPE)i * (DATA_TYPE)j) / (DATA_TYPE)N;
		B_[i * N + j] = ((DATA_TYPE)i * (DATA_TYPE)j + 1) / (DATA_TYPE)N;
	  }
	}
  }
};

class mvt : public BenchKernel {
 public:
  RT_NAME("mvt")
  explicit mvt(EDCL *edcl, int n = 4096) : BenchKernel(edcl) {
	N = n;
	a = edcl->createBuffer(N * N * sizeof(DATA_TYPE));
	x1 = edcl->createBuffer(N * sizeof(DATA_TYPE));
	x2 = edcl->createBuffer(N * sizeof(DATA_TYPE));
	y1 = edcl->createBuffer(N * sizeof(DATA_TYPE));
	y2 = edcl->createBuffer(N * sizeof(DATA_TYPE));
	init_arrays(PTR(a), PTR(x1), PTR(x2), PTR(y1), PTR(y2));
	kernelSource =
#include "Kernels/mvt.cl"
	program_ = edcl->createProgram(&kernelSource);
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = 1;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(2);
	EDKernel k1 = edcl_->createKernel(program_, "mvt_kernel1", 4);
	k1->configKernel(1, globalWorkSize, localWorkSize, a, x1, y1, &N);
	EDKernel k2 = edcl_->createKernel(program_, "mvt_kernel2", 4);
	k2->configKernel(1, globalWorkSize, localWorkSize, a, x2, y2, &N);
	aks->addKernel(k1, k2);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *a_ = PTR(a);
	DATA_TYPE *x1_ = PTR(x1);
	DATA_TYPE *x2_ = PTR(x2);
	DATA_TYPE *y1_ = PTR(y1);
	DATA_TYPE *y2_ = PTR(y2);
	int i, j, k, l;
	double timer = getCurrentTime();
#pragma omp parallel shared(a_, x1_, x2_, y1_, y2_, N) private(i, j, k, l)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		  x1_[i] = x1_[i] + a_[i * N + j] * y1_[j];
		}
	  }
#pragma omp for schedule(static)
	  for (k = 0; k < N; k++) {
		for (l = 0; l < N; l++) {
		  x2_[k] = x2_[k] + a_[k * N + l] * y2_[l];
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  EDBuffer a, x1, x2, y1, y2;
  int N;
  size_t DIM_LOCAL_WORK_GROUP_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_Y = 1;
  void init_arrays(DATA_TYPE *a_, DATA_TYPE *x1_, DATA_TYPE *x2_, DATA_TYPE *y_1, DATA_TYPE *y_2) const {
	int i, j;
	for (i = 0; i < N; i++) {
	  x1_[i] = 0.0;
	  x2_[i] = 0.0;
	  y_1[i] = 0.0;
	  y_2[i] = 0.0;
	  for (j = 0; j < N; j++) {
		a_[i * N + j] = (DATA_TYPE)(i + j + 1.0) / (DATA_TYPE)N;
	  }
	}
  }
};

class GramSchmidt : public BenchKernel {
 public:
  RT_NAME("GramSchmidt")
  explicit GramSchmidt(EDCL *edcl, int m = 1024, int n = 1024) : BenchKernel(edcl) {
	M = m;
	N = n;
	A = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	R = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	Q = edcl->createBuffer(M * N * sizeof(DATA_TYPE));
	init_array(PTR(A));
	kernelSource =
#include "Kernels/gramschmidt.cl"
	program_ = edcl->createProgram(&kernelSource);
	localWorkSize[0] = DIM_THREAD_BLOCK_X;
	localWorkSize[1] = DIM_THREAD_BLOCK_Y;
	globalWorkSizeKernel1[0] = DIM_THREAD_BLOCK_X;
	globalWorkSizeKernel1[1] = DIM_THREAD_BLOCK_Y;
	globalWorkSizeKernel2[0] = (size_t)ceil(((float)N) / ((float)DIM_THREAD_BLOCK_X)) * DIM_THREAD_BLOCK_X;
	globalWorkSizeKernel2[1] = 1;
	globalWorkSizeKernel3[0] = (size_t)ceil(((float)N) / ((float)DIM_THREAD_BLOCK_X)) * DIM_THREAD_BLOCK_X;
	globalWorkSizeKernel3[1] = 1;
	K = static_cast<int *>(calloc(N, sizeof(*K)));
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(3 * N);
	for (int k = 0; k < N; ++k) {
	  K[k] = k;
	  EDKernel k1 = edcl_->createKernel(program_, "gramschmidt_kernel1", 6);
	  k1->configKernel(1, globalWorkSizeKernel1, localWorkSize, A, R, Q, &K[k], &M, &N);
	  EDKernel k2 = edcl_->createKernel(program_, "gramschmidt_kernel2", 6);
	  k2->configKernel(1, globalWorkSizeKernel2, localWorkSize, A, R, Q, &K[k], &M, &N);
	  EDKernel k3 = edcl_->createKernel(program_, "gramschmidt_kernel3", 6);
	  k3->configKernel(1, globalWorkSizeKernel3, localWorkSize, A, R, Q, &K[k], &M, &N);
	  aks->addKernel(k1, k2, k3);
	}
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *R_ = PTR(R);
	DATA_TYPE *Q_ = PTR(Q);
	int i, j, k;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, R_, Q_, M, N) private(i, j, k)
	{
	  DATA_TYPE nrm;
#pragma omp for schedule(static)
	  for (k = 0; k < N; k++) {
		nrm = 0;
		for (i = 0; i < M; i++) {
		  nrm += A_[i * N + k] * A_[i * N + k];
		}

		R_[k * N + k] = sqrt(nrm);
		for (i = 0; i < M; i++) {
		  Q_[i * N + k] = A_[i * N + k] / R_[k * N + k];
		}

		for (j = k + 1; j < N; j++) {
		  R_[k * N + j] = 0;
		  for (i = 0; i < M; i++) {
			R_[k * N + j] += Q_[i * N + k] * A_[i * N + j];
		  }
		  for (i = 0; i < M; i++) {
			A_[i * N + j] = A_[i * N + j] - Q_[i * N + k] * R_[k * N + j];
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  int M, N;
  int *K;
  EDBuffer A, R, Q;
  size_t DIM_THREAD_BLOCK_X = 256;
  size_t DIM_THREAD_BLOCK_Y = 1;
  size_t globalWorkSizeKernel1[2]{};
  size_t globalWorkSizeKernel2[2]{};
  size_t globalWorkSizeKernel3[2]{};
  void init_array(DATA_TYPE *A_) const {
	int i, j;
	for (i = 0; i < M; i++) {
	  for (j = 0; j < N; j++) {
		A_[i * N + j] = ((DATA_TYPE)(i + 1) * (DATA_TYPE)(j + 1)) / (DATA_TYPE)(M + 1);
	  }
	}
  }
};

class Gesummv : public BenchKernel {
 public:
  RT_NAME("gesummv")
  explicit Gesummv(EDCL *edcl, int n = 4096) : BenchKernel(edcl) {
	N = n;
	A = edcl->createBuffer(N * N * sizeof(DATA_TYPE));
	B = edcl->createBuffer(N * N * sizeof(DATA_TYPE));
	x = edcl->createBuffer(N * sizeof(DATA_TYPE));
	y = edcl->createBuffer(N * sizeof(DATA_TYPE));
	tmp = edcl->createBuffer(N * sizeof(DATA_TYPE));
	init(PTR(A), PTR(x));
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
	workDim = 1;
	numArgs = 8;
	kernelSource =
#include "Kernels/gesummv.cl"
	program_ = edcl->createProgram(&kernelSource);
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(1);
	EDKernel kernel = edcl_->createKernel(program_, "gesummv_kernel", numArgs);
	kernel->configKernel(workDim, globalWorkSize, localWorkSize,
						 A, B, x, y, tmp, &alpha, &beta, &N);
	aks->addKernel(kernel);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	DATA_TYPE *x_ = PTR(x);
	DATA_TYPE *y_ = PTR(y);
	DATA_TYPE *tmp_ = PTR(tmp);
	int i, j;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, x_, y_, tmp_, N) private(i, j)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < N; i++) {
		tmp_[i] = 0;
		y_[i] = 0;
		for (j = 0; j < N; j++) {
		  tmp_[i] = A_[i * N + j] * x_[j] + tmp_[i];
		  y_[i] = B_[i * N + j] * x_[j] + y_[i];
		}

		y_[i] = (DATA_TYPE)alpha * tmp_[i] + (DATA_TYPE)beta * y_[i];
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  EDBuffer A, B, x, y, tmp;
  int N;
  int alpha = 1;
  int beta = 1;
  size_t DIM_LOCAL_WORK_GROUP_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_Y = 1;
  void init(DATA_TYPE *A_, DATA_TYPE *x_) const {
	int i, j;
	for (i = 0; i < N; i++) {
	  x_[i] = ((DATA_TYPE)i) / (DATA_TYPE)N;
	  for (j = 0; j < N; j++) {
		A_[i * N + j] = (DATA_TYPE)(i * j) / (DATA_TYPE)N;
	  }
	}
  }
};

class GEMM : public BenchKernel {
 public:
  RT_NAME("gemm")
  explicit GEMM(EDCL *edcl, int ni = 512, int nj = 512, int nk = 512) : BenchKernel(edcl) {
	NI = ni;
	NJ = nj;
	NK = nk;
	A = edcl->createBuffer(NI * NK * sizeof(DATA_TYPE));
	B = edcl->createBuffer(NK * NJ * sizeof(DATA_TYPE));
	C = edcl->createBuffer(NI * NJ * sizeof(DATA_TYPE));
	init(PTR(A), PTR(B), PTR(C));
	kernelSource =
#include "Kernels/gemm.cl"
	program_ = edcl->createProgram(&kernelSource);
	workDim = 2;
	numArgs = 8;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NJ) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)NI) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	DATA_TYPE *C_ = PTR(C);
	int i, j, k;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, C_, NI, NJ, NK) private(i, j, k)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
		  C_[i * NJ + j] *= (DATA_TYPE)beta;
		  for (k = 0; k < NK; ++k) {
			C_[i * NJ + j] += (DATA_TYPE)alpha * A_[i * NK + k] * B_[k * NJ + j];
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(1);
	EDKernel gemm = edcl_->createKernel(program_, "gemm", numArgs);
	gemm->configKernel(workDim, globalWorkSize, localWorkSize, A, B, C, &alpha, &beta, &NI, &NJ, &NK);
	aks->addKernel(gemm);
	return aks;
  }

 private:
  EDBuffer A, B, C;
  int NI, NJ, NK;
  int alpha = 32412;
  int beta = 2123;
  size_t DIM_LOCAL_WORK_GROUP_X = 32;
  size_t DIM_LOCAL_WORK_GROUP_Y = 8;
  void init(DATA_TYPE *A_, DATA_TYPE *B_, DATA_TYPE *C_) const {
	int i, j;
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NK; j++) {
		A_[i * NK + j] = (DATA_TYPE)(i * j) / (DATA_TYPE)NI;
	  }
	}
	for (i = 0; i < NK; i++) {
	  for (j = 0; j < NJ; j++) {
		B_[i * NJ + j] = (DATA_TYPE)(i * j + 1) / (DATA_TYPE)NJ;
	  }
	}
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NJ; j++) {
		C_[i * NJ + j] = ((DATA_TYPE)i * (DATA_TYPE)j + 2) / (DATA_TYPE)NJ;
	  }
	}
  }

};

class FDTD2D : public BenchKernel {
 public:
  RT_NAME("FDTD2D")
  explicit FDTD2D(EDCL *edcl, int tmax = 500, int nx = 2048, int ny = 2048) : BenchKernel(edcl) {
	TMAX = tmax;
	NX = nx;
	NY = ny;
	_fict_ = edcl->createBuffer(TMAX * sizeof(DATA_TYPE));
	ex = edcl->createBuffer(NX * (NY + 1) * sizeof(DATA_TYPE));
	ey = edcl->createBuffer((NX + 1) * NY * sizeof(DATA_TYPE));
	hz = edcl->createBuffer(NX * NY * sizeof(DATA_TYPE));
	T = static_cast<int *>(calloc(TMAX, sizeof(*T)));
	init_arrays(PTR(_fict_), PTR(ex), PTR(ey), PTR(hz));
	kernelSource =
#include "Kernels/fdtd2d.cl"
	program_ = edcl->createProgram(&kernelSource);
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NY) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)NX) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
	workDim = 2;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(TMAX * 3);
	for (int t = 0; t < TMAX; ++t) {
	  T[t] = t;
	  EDKernel k1 = edcl_->createKernel(program_, "fdtd_kernel1", 7);
	  k1->configKernel(workDim, globalWorkSize, localWorkSize, _fict_, ex, ey, hz, &T[t], &NX, &NY);
	  EDKernel k2 = edcl_->createKernel(program_, "fdtd_kernel2", 5);
	  k2->configKernel(workDim, globalWorkSize, localWorkSize, ex, ey, hz, &NX, &NY);
	  EDKernel k3 = edcl_->createKernel(program_, "fdtd_kernel3", 5);
	  k3->configKernel(workDim, globalWorkSize, localWorkSize, ex, ey, hz, &NX, &NY);
	  aks->addKernel(k1, k2, k3);
	}
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *fict_ = PTR(_fict_);
	DATA_TYPE *ex_ = PTR(ex);
	DATA_TYPE *ey_ = PTR(ey);
	DATA_TYPE *hz_ = PTR(hz);
	int t, i, j;
	double timer = getCurrentTime();
#pragma omp parallel shared(fict_, ex_, ey_, hz_, TMAX, NY, NX) private(i, j, t)
	{
#pragma omp for schedule(static)
	  for (t = 0; t < TMAX; t++) {
		for (j = 0; j < NY; j++) {
		  ey_[0 * NY + j] = fict_[t];
		}
		for (i = 1; i < NX; i++) {
		  for (j = 0; j < NY; j++) {
			ey_[i * NY + j] = ey_[i * NY + j] - 0.5 * (hz_[i * NY + j] - hz_[(i - 1) * NY + j]);
		  }
		}
		for (i = 0; i < NX; i++) {
		  for (j = 1; j < NY; j++) {
			ex_[i * (NY + 1) + j] = ex_[i * (NY + 1) + j] - 0.5 * (hz_[i * NY + j] - hz_[i * NY + (j - 1)]);
		  }
		}
		for (i = 0; i < NX; i++) {
		  for (j = 0; j < NY; j++) {
			hz_[i * NY + j] = hz_[i * NY + j]
				- 0.7 * (ex_[i * (NY + 1) + (j + 1)] - ex_[i * (NY + 1) + j] + ey_[(i + 1) * NY + j] - ey_[i * NY + j]);
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

  ~FDTD2D() {
	delete T;
  }
 private:
  EDBuffer _fict_;
  EDBuffer ex;
  EDBuffer ey;
  EDBuffer hz;
  int *T;
  int TMAX, NX, NY;
  size_t DIM_LOCAL_WORK_GROUP_X = 32;
  size_t DIM_LOCAL_WORK_GROUP_Y = 8;
  void init_arrays(DATA_TYPE *fict, DATA_TYPE *ex_, DATA_TYPE *ey_, DATA_TYPE *hz_) const {
	int i, j;
	for (i = 0; i < TMAX; i++) {
	  fict[i] = (DATA_TYPE)i;
	}
	for (i = 0; i < NX; i++) {
	  for (j = 0; j < NY; j++) {
		ex_[i * NY + j] = ((DATA_TYPE)i * (DATA_TYPE)(j + 1) + 1) / (DATA_TYPE)NX;
		ey_[i * NY + j] = ((DATA_TYPE)(i - 1) * (DATA_TYPE)(j + 2) + 2) / (DATA_TYPE)NX;
		hz_[i * NY + j] = ((DATA_TYPE)(i - 9) * (DATA_TYPE)(j + 4) + 3) / (DATA_TYPE)NX;
	  }
	}
  }
};

#define FLOAT_N 3214212.01
#define EPS 0.005
class Covariance : public BenchKernel {
 public:
  RT_NAME("Covariance")
  explicit Covariance(EDCL *edcl, int m = 1024, int n = 1024) : BenchKernel(edcl) {
	M = m;
	N = n;
	data = edcl->createBuffer((M + 1) * (N + 1) * sizeof(DATA_TYPE));
	symmat = edcl->createBuffer((M + 1) * (N + 1) * sizeof(DATA_TYPE));
	mean = edcl->createBuffer((M + 1) * sizeof(DATA_TYPE));
	init_arrays(PTR(data));
	kernelSource =
#include "Kernels/covariance.cl"
	program_ = edcl->createProgram(&kernelSource);
	localWorkSize_Kernel1[0] = DIM_LOCAL_WORK_GROUP_KERNEL_1_X;
	localWorkSize_Kernel1[1] = DIM_LOCAL_WORK_GROUP_KERNEL_1_Y;
	globalWorkSize_Kernel1[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_1_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_1_X;
	globalWorkSize_Kernel1[1] = 1;
	localWorkSize_Kernel2[0] = DIM_LOCAL_WORK_GROUP_KERNEL_2_X;
	localWorkSize_Kernel2[1] = DIM_LOCAL_WORK_GROUP_KERNEL_2_Y;
	globalWorkSize_Kernel2[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_2_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_2_X;
	globalWorkSize_Kernel2[1] =
		(size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_2_Y)) * DIM_LOCAL_WORK_GROUP_KERNEL_2_Y;
	localWorkSize_Kernel3[0] = DIM_LOCAL_WORK_GROUP_KERNEL_3_X;
	localWorkSize_Kernel3[1] = DIM_LOCAL_WORK_GROUP_KERNEL_3_Y;
	globalWorkSize_Kernel3[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_3_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_3_X;
	globalWorkSize_Kernel3[1] = 1;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(3);
	EDKernel kernel_mean = edcl_->createKernel(program_, k1Name, 5);
	kernel_mean->configKernel(1, globalWorkSize_Kernel1, localWorkSize_Kernel1, mean, data, &float_n, &M, &N);
	EDKernel kernel_reduce = edcl_->createKernel(program_, k2Name, 4);
	kernel_reduce->configKernel(2, globalWorkSize_Kernel2, localWorkSize_Kernel2, mean, data, &M, &N);
	EDKernel kernel_covar = edcl_->createKernel(program_, k3Name, 4);
	kernel_covar->configKernel(1, globalWorkSize_Kernel3, localWorkSize_Kernel3, symmat, data, &M, &N);
	aks->addKernel(kernel_mean, kernel_reduce, kernel_covar);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *mean_ = PTR(mean);
	DATA_TYPE *symmat_ = PTR(symmat);
	DATA_TYPE *data_ = PTR(data);
	int i, j, j1, j2;
	double timer = getCurrentTime();
#pragma omp parallel shared(mean, symmat_, data, M, N) private(i, j, j1, j2)
	{
	  /* Determine mean of column vectors of input data matrix */
#pragma omp for schedule(static)
	  for (j = 1; j <= M; j++) {
		mean_[j] = 0.0;
		for (i = 1; i <= N; i++) {
		  mean_[j] += data_[i * (M + 1) + j];
		}
		mean_[j] /= float_n;
	  }

	  /* Center the column vectors. */
#pragma omp for schedule(static)
	  for (i = 1; i <= N; i++) {
		for (j = 1; j <= M; j++) {
		  data_[i * (M + 1) + j] -= mean_[j];
		}
	  }

	  /* Calculate the m * m covariance matrix. */
#pragma omp for schedule(static)
	  for (j1 = 1; j1 <= M; j1++) {
		for (j2 = j1; j2 <= M; j2++) {
		  symmat_[j1 * (M + 1) + j2] = 0.0;
		  for (i = 1; i <= N; i++) {
			symmat_[j1 * (M + 1) + j2] += data_[i * (M + 1) + j1] * data_[i * (M + 1) + j2];
		  }
		  symmat_[j2 * (M + 1) + j1] = symmat_[j1 * (M + 1) + j2];
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  int M, N;
  EDBuffer data, symmat, mean;
  DATA_TYPE float_n = FLOAT_N;
  std::string k1Name = "mean_kernel";
  std::string k2Name = "reduce_kernel";
  std::string k3Name = "covar_kernel";
  /* Thread block dimensions for kernel 1*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_1_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_1_Y = 1;

/* Thread block dimensions for kernel 2*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_2_X = 32;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_2_Y = 8;

/* Thread block dimensions for kernel 3*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_3_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_3_Y = 1;

  size_t localWorkSize_Kernel1[2]{};
  size_t localWorkSize_Kernel2[2]{};
  size_t localWorkSize_Kernel3[2]{};
  size_t globalWorkSize_Kernel1[2]{};
  size_t globalWorkSize_Kernel2[2]{};
  size_t globalWorkSize_Kernel3[2]{};

  void init_arrays(DATA_TYPE *data_) const {
	int i, j;
	for (i = 0; i < M; i++) {
	  for (j = 0; j < N; j++) {
		data_[i * (N + 1) + j] = (DATA_TYPE)(i * j) / (DATA_TYPE)M;
	  }
	}
  }
};

class Correlation : public BenchKernel {
 public:
  RT_NAME("Correlation")
  explicit Correlation(EDCL *edcl, int m = 1024, int n = 1024) : BenchKernel(edcl) {
	M = m;
	N = n;
	data = edcl->createBuffer((M + 1) * (N + 1) * sizeof(DATA_TYPE));
	mean = edcl->createBuffer((M + 1) * sizeof(DATA_TYPE));
	stddev = edcl->createBuffer((M + 1) * sizeof(DATA_TYPE));
	symmat = edcl->createBuffer((M + 1) * (N + 1) * sizeof(DATA_TYPE));
	init_arrays(PTR(data));
	kernelSource =
#include "Kernels/correlation.cl"
	program_ = edcl_->createProgram(&kernelSource);
	workDim = 2;
	localWorkSize_Kernel1[0] = DIM_LOCAL_WORK_GROUP_KERNEL_1_X;
	localWorkSize_Kernel1[1] = DIM_LOCAL_WORK_GROUP_KERNEL_1_Y;
	globalWorkSize_Kernel1[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_1_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_1_X;
	globalWorkSize_Kernel1[1] = 1;

	localWorkSize_Kernel2[0] = DIM_LOCAL_WORK_GROUP_KERNEL_2_X;
	localWorkSize_Kernel2[1] = DIM_LOCAL_WORK_GROUP_KERNEL_2_Y;
	globalWorkSize_Kernel2[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_2_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_2_X;
	globalWorkSize_Kernel2[1] = 1;

	localWorkSize_Kernel3[0] = DIM_LOCAL_WORK_GROUP_KERNEL_3_X;
	localWorkSize_Kernel3[1] = DIM_LOCAL_WORK_GROUP_KERNEL_3_Y;
	globalWorkSize_Kernel3[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_3_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_3_X;
	globalWorkSize_Kernel3[1] =
		(size_t)ceil(((float)N) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_3_Y)) * DIM_LOCAL_WORK_GROUP_KERNEL_3_Y;

	localWorkSize_Kernel4[0] = DIM_LOCAL_WORK_GROUP_KERNEL_4_X;
	localWorkSize_Kernel4[1] = DIM_LOCAL_WORK_GROUP_KERNEL_4_Y;
	globalWorkSize_Kernel4[0] =
		(size_t)ceil(((float)M) / ((float)DIM_LOCAL_WORK_GROUP_KERNEL_4_X)) * DIM_LOCAL_WORK_GROUP_KERNEL_4_X;
	globalWorkSize_Kernel4[1] = 1;

  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(4);
	EDKernel kernel_mean = edcl_->createKernel(program_, "mean_kernel", numArgs_mean);
	EDKernel kernel_std = edcl_->createKernel(program_, "std_kernel", numArgs_std);
	EDKernel kernel_reduce = edcl_->createKernel(program_, "reduce_kernel", numArgs_reduce);
	EDKernel kernel_corr = edcl_->createKernel(program_, "corr_kernel", numArgs_corr);
	kernel_mean->configKernel(1, globalWorkSize_Kernel1, localWorkSize_Kernel1, mean, data, &float_n, &M, &N);
	kernel_std->configKernel(1, globalWorkSize_Kernel2, localWorkSize_Kernel2,
							 mean, stddev, data, &float_n, &eps, &M, &N);
	kernel_reduce->configKernel(2, globalWorkSize_Kernel3, localWorkSize_Kernel3, mean, stddev, data, &float_n, &M, &N);
	kernel_corr->configKernel(1, globalWorkSize_Kernel4, localWorkSize_Kernel4, symmat, data, &M, &N);
	aks->addKernel(kernel_mean, kernel_std, kernel_reduce, kernel_corr);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *data_ = PTR(data);
	DATA_TYPE *mean_ = PTR(mean);
	DATA_TYPE *stddev_ = PTR(stddev);
	DATA_TYPE *symmat_ = PTR(symmat);
	int i, j, j1, j2;
	double timer = getCurrentTime();
#pragma omp parallel shared(data_, mean_, stddev_, symmat_, M, N) private(i, j, j1, j2)
	{
#pragma omp for schedule(static)
	  // Determine mean of column vectors of input data matrix
	  for (j = 1; j <= M; j++) {
		mean_[j] = 0.0;
		for (i = 1; i <= N; i++) {
		  mean_[j] += data_[i * (M + 1) + j];
		}
		mean_[j] /= (DATA_TYPE)FLOAT_N;
	  }

	  // Determine standard deviations of column vectors of data matrix.
#pragma omp for schedule(static)
	  for (j = 1; j <= M; j++) {
		stddev_[j] = 0.0;
		for (i = 1; i <= N; i++) {
		  stddev_[j] += (data_[i * (M + 1) + j] - mean_[j]) * (data_[i * (M + 1) + j] - mean_[j]);
		}
		stddev_[j] /= FLOAT_N;
		stddev_[j] = sqrt_of_array_cell(stddev_, j);
		stddev_[j] = stddev_[j] <= EPS ? 1.0 : stddev_[j];
	  }

	  // Center and reduce the column vectors.
#pragma omp for schedule(static)
	  for (i = 1; i <= N; i++) {
		for (j = 1; j <= M; j++) {
		  data_[i * (M + 1) + j] -= mean_[j];
		  data_[i * (M + 1) + j] /= sqrt(FLOAT_N);
		  data_[i * (M + 1) + j] /= stddev_[j];
		}
	  }

	  // Calculate the m * m correlation matrix.
#pragma omp for schedule(static)
	  for (j1 = 1; j1 <= M - 1; j1++) {
		symmat_[j1 * (M + 1) + j1] = 1.0;

		for (j2 = j1 + 1; j2 <= M; j2++) {
		  symmat_[j1 * (M + 1) + j2] = 0.0;

		  for (i = 1; i <= N; i++) {
			symmat_[j1 * (M + 1) + j2] += (data_[i * (M + 1) + j1] * data_[i * (M + 1) + j2]);
		  }

		  symmat_[j2 * (M + 1) + j1] = symmat_[j1 * (M + 1) + j2];
		}
	  }

	  symmat_[M * (M + 1) + M] = 1.0;
	}
	return getCurrentTime() - timer;
  }

 private:
  DATA_TYPE float_n = FLOAT_N;
  DATA_TYPE eps = EPS;
  int M, N;
  EDBuffer data, mean, stddev, symmat;
  /* Thread block dimensions for kernel 1*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_1_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_1_Y = 1;

/* Thread block dimensions for kernel 2*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_2_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_2_Y = 1;

/* Thread block dimensions for kernel 3*/
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_3_X = 32;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_3_Y = 8;

/* Thread block dimensions for kernel 4*/
//  size_t DIM_LOCAL_WORK_GROUP_KERNEL_4_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_4_X = 256;
  size_t DIM_LOCAL_WORK_GROUP_KERNEL_4_Y = 1;
  size_t localWorkSize_Kernel1[2]{};
  size_t localWorkSize_Kernel2[2]{};
  size_t localWorkSize_Kernel3[2]{};
  size_t localWorkSize_Kernel4[2]{};
  size_t globalWorkSize_Kernel1[2]{};
  size_t globalWorkSize_Kernel2[2]{};
  size_t globalWorkSize_Kernel3[2]{};
  size_t globalWorkSize_Kernel4[2]{};
  int numArgs_mean = 5;
  int numArgs_std = 7;
  int numArgs_reduce = 6;
  int numArgs_corr = 4;
  void init_arrays(DATA_TYPE *data_) const {
	int i, j;
	for (i = 0; i <= M; i++) {
	  for (j = 0; j <= N; j++) {
		data_[i * N + j] = ((DATA_TYPE)i * (DATA_TYPE)j) / ((DATA_TYPE)M + 1);
	  }
	}
  }
};

/* Thread block dimensions */
#define DIM_LOCAL_WORK_GROUP_X 256
#define DIM_LOCAL_WORK_GROUP_Y 1
#ifndef M_PI
#define M_PI 3.14159265358979323846	/* pi */
#endif

class BICG : public BenchKernel {
 public:
  RT_NAME("bicg")
  explicit BICG(EDCL *edcl, int nx = 4096, int ny = 4096) : BenchKernel(edcl) {
	NX = nx;
	NY = ny;
	a = edcl->createBuffer(sizeof(DATA_TYPE) * NX * NY);
	r = edcl->createBuffer(sizeof(DATA_TYPE) * NX);
	s = edcl->createBuffer(sizeof(DATA_TYPE) * NX);
	p = edcl->createBuffer(sizeof(DATA_TYPE) * NX);
	q = edcl->createBuffer(sizeof(DATA_TYPE) * NX);
	init_array(PTR(a), PTR(p), PTR(r));
	kernelSource =
#include "Kernels/bicg.cl"
	program_ = edcl->createProgram(&kernelSource);
	numArgs = 5;
	workDim = 1;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NX) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = 1;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(2);
	EDKernel kernel1 = edcl_->createKernel(program_, kernel1Name, numArgs);
	kernel1->configKernel(workDim, globalWorkSize, localWorkSize,
						  a, p, q, &NX, &NY);
	EDKernel kernel2 = edcl_->createKernel(program_, kernel2Name, numArgs);
	kernel2->configKernel(workDim, globalWorkSize, localWorkSize,
						  a, r, s, &NX, &NY);
	aks->addKernel(kernel1, kernel2);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	DATA_TYPE *A_ = PTR(a);
	DATA_TYPE *r_ = PTR(r);
	DATA_TYPE *s_ = PTR(s);
	DATA_TYPE *p_ = PTR(p);
	DATA_TYPE *q_ = PTR(q);
	int i, j;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, r_, s_, p_, q_, NX, NY) private(i, j)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < NY; i++) {
		s_[i] = 0.0;
	  }
#pragma omp for schedule(static)
	  for (i = 0; i < NX; i++) {
		q_[i] = 0.0;
		for (j = 0; j < NY; j++) {
		  s_[j] = s_[j] + r_[i] * A_[i * NY + j];
		  q_[i] = q_[i] + A_[i * NY + j] * p_[j];
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  std::string kernel1Name = "bicgKernel1";
  std::string kernel2Name = "bicgKernel2";
  EDBuffer a, r, s, p, q;
  int NX, NY;
  void init_array(DATA_TYPE *A_, DATA_TYPE *p_, DATA_TYPE *r_) const {
	int i, j;
	for (i = 0; i < NX; i++) {
	  r_[i] = i * M_PI;
	  for (j = 0; j < NY; j++) {
		A_[i * NY + j] = ((DATA_TYPE)i * (DATA_TYPE)j) / (DATA_TYPE)NX;
	  }
	}
	for (i = 0; i < NY; i++) {
	  p_[i] = i * M_PI;
	}
  }
};

class ATAX : public BenchKernel {
 public:
  RT_NAME("atax")
  explicit ATAX(EDCL *edcl, int nx = 4096, int ny = 4096) : BenchKernel(edcl) {
	NX = nx;
	NY = ny;
	A = edcl->createBuffer(NX * NY * sizeof(DATA_TYPE));
	x = edcl->createBuffer(NX * sizeof(DATA_TYPE));
	y = edcl->createBuffer(NY * sizeof(DATA_TYPE));
	tmp = edcl->createBuffer(NX * sizeof(DATA_TYPE));
	init_array(PTR(x), PTR(A));
	kernelSource =
#include "Kernels/atax.cl"
	program_ = edcl->createProgram(&kernelSource);
	numArgs = 5;
	workDim = 1;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NX) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = 1;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(2);
	EDKernel kernel1 = edcl_->createKernel(program_, kernel1Name, numArgs);
	kernel1->configKernel(workDim, globalWorkSize, localWorkSize,
						  A, x, tmp, &NX, &NY);
	EDKernel kernel2 = edcl_->createKernel(program_, kernel2Name, numArgs);
	kernel2->configKernel(workDim, globalWorkSize, localWorkSize,
						  A, x, tmp, &NX, &NY);
	aks->addKernel(kernel1, kernel2);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *x_ = PTR(x);
	DATA_TYPE *y_ = PTR(y);
	DATA_TYPE *tmp_ = PTR(tmp);
	int i, j;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, x_, y_, tmp_, NX, NY) private(i, j)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < NY; i++) {
		y_[i] = 0;
	  }
#pragma omp for schedule(static)
	  for (i = 0; i < NX; i++) {
		tmp_[i] = 0;
		for (j = 0; j < NY; j++) {
		  tmp_[i] = tmp_[i] + A_[i * NY + j] * x_[j];
		}
		for (j = 0; j < NY; j++) {
		  y_[j] = y_[j] + A_[i * NY + j] * tmp_[i];
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  std::string kernel1Name = "atax_kernel1";
  std::string kernel2Name = "atax_kernel2";
  int NX, NY;
  EDBuffer A, x, y, tmp;
  void init_array(DATA_TYPE *x_, DATA_TYPE *A_) const {
	int i, j;
	for (i = 0; i < NX; i++) {
	  x_[i] = i * M_PI;
	  for (j = 0; j < NY; j++) {
		A_[i * NY + j] = ((DATA_TYPE)i * (DATA_TYPE)(j)) / (DATA_TYPE)NX;
	  }
	}
  }
};
#undef DIM_LOCAL_WORK_GROUP_X
#undef DIM_LOCAL_WORK_GROUP_Y

/* Thread block dimensions */
#define DIM_LOCAL_WORK_GROUP_X 32
#define DIM_LOCAL_WORK_GROUP_Y 8
class MM2 : public BenchKernel {
 public:
  RT_NAME("2mm")
  explicit MM2(EDCL *edcl, int ni = 1024, int nj = 1024, int nk = 1024, int nl = 1024) : BenchKernel(edcl) {
	NI = ni;
	NJ = nj;
	NK = nk;
	NL = nl;
	C = edcl->createBuffer(NI * NJ * sizeof(float));
	A = edcl->createBuffer(NI * NK * sizeof(float));
	B = edcl->createBuffer(NK * NJ * sizeof(float));
	D = edcl->createBuffer(NJ * NL * sizeof(float));
	E = edcl->createBuffer(NI * NL * sizeof(float));
	init_array(PTR(A), PTR(B), PTR(C), PTR(D));
	kernelSource =
#include "Kernels/2mm.cl"
	program_ = edcl_->createProgram(&kernelSource);
	numArgs = 6;
	workDim = 2;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NI) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)NL) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(2);
	EDKernel kernel1 = edcl_->createKernel(program_, kernel1Name, numArgs);
	kernel1->configKernel(workDim, globalWorkSize, localWorkSize,
						  A, B, C, &NI, &NK, &NJ);
	EDKernel kernel2 = edcl_->createKernel(program_, kernel2Name, numArgs);
	kernel2->configKernel(workDim, globalWorkSize, localWorkSize,
						  C, D, E, &NI, &NJ, &NL);
	aks->addKernel(kernel1, kernel2);
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	DATA_TYPE *C_ = PTR(C);
	DATA_TYPE *D_ = PTR(D);
	DATA_TYPE *E_ = PTR(E);
	int i, j, k;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, C_, D_, E_, NI, NJ, NK, NL) private(i, j, k)
	{
#pragma omp for schedule(static)
	  for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
		  for (k = 0; k < NK; ++k) {
			C_[i * NJ + j] += A_[i * NK + k] * B_[k * NJ + j];
		  }
		}
	  }
#pragma omp for schedule(static)
	  for (i = 0; i < NI; i++) {
		for (j = 0; j < NL; j++) {
		  for (k = 0; k < NJ; ++k) {
			E_[i * NL + j] += C_[i * NJ + k] * D_[k * NL + j];
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  std::string kernel1Name = "mm2_kernel1";
  std::string kernel2Name = "mm2_kernel2";
  EDBuffer A, B, C, D, E;
  int NI, NJ, NK, NL;
  void init_array(DATA_TYPE *A_, DATA_TYPE *B_, DATA_TYPE *C_, DATA_TYPE *D_) const {
	int i, j;

	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NK; j++) {
		A_[i * NI + j] = ((DATA_TYPE)i * (DATA_TYPE)j) / (DATA_TYPE)NI;
	  }
	}

	for (i = 0; i < NK; i++) {
	  for (j = 0; j < NJ; j++) {
		B_[i * NK + j] = ((DATA_TYPE)i * (DATA_TYPE)(j + 1)) / (DATA_TYPE)NJ;
	  }
	}

	for (i = 0; i < NL; i++) {
	  for (j = 0; j < NJ; j++) {
		C_[i * NL + j] = ((DATA_TYPE)i * (DATA_TYPE)(j + 3)) / (DATA_TYPE)NL;
	  }
	}

	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NL; j++) {
		D_[i * NL + j] = ((DATA_TYPE)i * (DATA_TYPE)(j + 2)) / (DATA_TYPE)NK;
	  }
	}
  }
};

class Conv3D : public BenchKernel {
 public:
  RT_NAME("Conv3D")
  explicit Conv3D(EDCL *edcl, int ni = 256, int nj = 256, int nk = 256) : BenchKernel(edcl) {
	kernelSource =
#include "Kernels/conv3d.cl"
	program_ = edcl_->createProgram(&kernelSource);
	kernelName = "Convolution3D_kernel";
	numArgs = 6;
	NI = ni;
	NJ = nj;
	NK = nk;
	i_s = (int *)calloc(NI, sizeof(*i_s));
	A = edcl->createBuffer(NI * NJ * NK * sizeof(float));
	// init(A)
	for (int i = 0; i < NI; ++i) {
	  for (int j = 0; j < NJ; ++j) {
		for (int k = 0; k < NK; ++k) {
		  FLOAT_PTR(A)[i * (NK * NJ) + j * NK + k] = (DATA_TYPE)(i % 12 + 2 * (j % 7) + 3 * (k % 13));
		}
	  }
	}
	B = edcl->createBuffer(NI * NJ * NK * sizeof(float));
	workDim = 2;
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NK) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)NJ) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
  }

  AtomicKernelSet createBenchAKS() override {
	AtomicKernelSet aks = edcl_->createAtomicKernelSet(NI);
	for (int i = 1; i < NI - 1; ++i) {
	  i_s[i] = i;
	  EDKernel kernel = edcl_->createKernel(program_, kernelName, numArgs);
	  kernel->configKernel(workDim, globalWorkSize, localWorkSize, A, B, &NI, &NJ, &NK, &i_s[i]);
	  aks->addKernel(kernel);
	}
	return aks;
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	int i, j, k;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, NI, NJ, NK) private(i, j, k)
	{
	  DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;
	  c11 = +2;
	  c21 = +5;
	  c31 = -8;
	  c12 = -3;
	  c22 = +6;
	  c32 = -9;
	  c13 = +4;
	  c23 = +7;
	  c33 = +10;
#pragma omp for schedule(static)
	  for (i = 1; i < NI - 1; ++i) {
		for (j = 1; j < NJ - 1; ++j) {
		  for (k = 1; k < NK - 1; ++k) {
			B_[i * (NK * NJ) + j * NK + k] = c11 * A_[(i - 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c13 * A_[(i + 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c21 * A_[(i - 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c23 * A_[(i + 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c31 * A_[(i - 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c33 * A_[(i + 1) * (NK * NJ) + (j - 1) * NK + (k - 1)]
				+ c12 * A_[(i + 0) * (NK * NJ) + (j - 1) * NK + (k + 0)]
				+ c22 * A_[(i + 0) * (NK * NJ) + (j + 0) * NK + (k + 0)]
				+ c32 * A_[(i + 0) * (NK * NJ) + (j + 1) * NK + (k + 0)]
				+ c11 * A_[(i - 1) * (NK * NJ) + (j - 1) * NK + (k + 1)]
				+ c13 * A_[(i + 1) * (NK * NJ) + (j - 1) * NK + (k + 1)]
				+ c21 * A_[(i - 1) * (NK * NJ) + (j + 0) * NK + (k + 1)]
				+ c23 * A_[(i + 1) * (NK * NJ) + (j + 0) * NK + (k + 1)]
				+ c31 * A_[(i - 1) * (NK * NJ) + (j + 1) * NK + (k + 1)]
				+ c33 * A_[(i + 1) * (NK * NJ) + (j + 1) * NK + (k + 1)];
		  }
		}
	  }
	}
	return getCurrentTime() - timer;
  }

  ~Conv3D() {
	delete i_s;
  }
 private:
  // Kernel args
  EDBuffer A, B;
  int NI, NJ, NK;
  int *i_s;
};

class Conv2D : public BenchKernel {
 public:
  RT_NAME("Conv2D")
  explicit Conv2D(EDCL *edcl, int ni = 4096, int nj = 4096) : BenchKernel(edcl) {
	kernelSource =
#include "Kernels/conv2d.cl"
	program_ = edcl_->createProgram(&kernelSource);
	kernelName = "Convolution2D_kernel";
	numArgs = 4;
	workDim = 2;
	NI = ni;
	NJ = nj;
	A = edcl_->createBuffer(NI * NJ * sizeof(float));
	// init(A)
	for (int i = 0; i < NI; ++i) {
	  for (int j = 0; j < NJ; ++j)
		FLOAT_PTR(A)[i * NJ + j] = randNumGenerator(0, 1);
	}
	B = edcl_->createBuffer(NI * NJ * sizeof(float));
	localWorkSize[0] = DIM_LOCAL_WORK_GROUP_X;
	localWorkSize[1] = DIM_LOCAL_WORK_GROUP_Y;
	globalWorkSize[0] = (size_t)ceil(((float)NI) / ((float)DIM_LOCAL_WORK_GROUP_X)) * DIM_LOCAL_WORK_GROUP_X;
	globalWorkSize[1] = (size_t)ceil(((float)NJ) / ((float)DIM_LOCAL_WORK_GROUP_Y)) * DIM_LOCAL_WORK_GROUP_Y;
  }

  AtomicKernelSet createBenchAKS() override {
	EDKernel kernel = edcl_->createKernel(program_, kernelName, numArgs);
	kernel->configKernel(workDim, globalWorkSize, localWorkSize, A, B, &NI, &NJ);
	return edcl_->wrapEDKernelInAtomicKernelSet(kernel);
  }

  double nativeKernel(int NumThreads) override {
	omp_set_num_threads(NumThreads);
	DATA_TYPE *A_ = PTR(A);
	DATA_TYPE *B_ = PTR(B);
	int i, j;
	double timer = getCurrentTime();
#pragma omp parallel shared(A_, B_, NI, NJ) private(i, j)
	{
	  DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;
	  c11 = +0.2;
	  c21 = +0.5;
	  c31 = -0.8;
	  c12 = -0.3;
	  c22 = +0.6;
	  c32 = -0.9;
	  c13 = +0.4;
	  c23 = +0.7;
	  c33 = +0.10;
#pragma omp for schedule(static)
	  for (i = 1; i < NI - 1; ++i) {
		for (j = 1; j < NJ - 1; ++j) {
		  B_[i * NJ + j] =
			  c11 * A_[(i - 1) * NJ + (j - 1)] + c12 * A_[(i + 0) * NJ + (j - 1)] + c13 * A_[(i + 1) * NJ + (j - 1)]
				  + c21 * A_[(i - 1) * NJ + (j + 0)] + c22 * A_[(i + 0) * NJ + (j + 0)]
				  + c23 * A_[(i + 1) * NJ + (j + 0)]
				  + c31 * A_[(i - 1) * NJ + (j + 1)] + c32 * A_[(i + 0) * NJ + (j + 1)]
				  + c33 * A_[(i + 1) * NJ + (j + 1)];
		}
	  }
	}
	return getCurrentTime() - timer;
  }

 private:
  EDBuffer A;
  EDBuffer B;
  int NI;
  int NJ;
};

#endif //EDCL_POLYBENCH_GPU_POLYBENCH_H_
