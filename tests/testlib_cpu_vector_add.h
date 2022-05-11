#ifndef EDCL_TESTS_TESTLIB_CPU_VECTOR_ADD_H_
#define EDCL_TESTS_TESTLIB_CPU_VECTOR_ADD_H_

void vecAddOptimized(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddNoneOptimized(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP8Opt(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP4Opt(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP8(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP4(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddOMP1(const float *va, const float *vb, float *vc, unsigned int n);

void vecAddNoneOptimized(const float *va,
						 const float *vb,
						 float *vc,
						 unsigned int n)__attribute__((optnone)) {
  int i;
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOptimized(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP8Opt(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(8)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP4Opt(const float *va, const float *vb, float *vc, unsigned int n) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(4)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP8(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(8)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP4(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(4)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

void vecAddOMP1(const float *va, const float *vb, float *vc, unsigned int n)__attribute__((optnone)) {
  int i;
#pragma omp parallel for default(none) private(i) shared(va, vb, vc, n) num_threads(1)
  for (i = 0; i < n; i++) {
	vc[i] = va[i] + vb[i];
  }
}

#endif //EDCL_TESTS_TESTLIB_CPU_VECTOR_ADD_H_
