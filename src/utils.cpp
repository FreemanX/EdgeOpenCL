#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#include <tar.h>
#include "utils.h"

unsigned int pow2(unsigned int v) {
  return static_cast<unsigned int>(1 << v);
}

__unused double sqr(double v) {
  return v * v;
}

__unused int compare_doubles(const void *a, const void *b) {
  const auto *v1 = (const double *)a;
  const auto *v2 = (const double *)b;
  if (*v1 == *v2) return 0;
  return *v1 > *v2 ? 1 : -1;
}

static pthread_mutex_t print_mx;
static bool print_mx_inited = false;
static void initDebugMX() {
  pthread_mutex_init(&print_mx, nullptr);
  print_mx_inited = true;
}

void flushed_printf(const char *format, ...) {
  if (!print_mx_inited) initDebugMX();
  pthread_mutex_lock(&print_mx);
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
  pthread_mutex_unlock(&print_mx);
}

void err_printf(const char *msgHead, int line, const char *what, ...) {
  if (!print_mx_inited) initDebugMX();
  pthread_mutex_lock(&print_mx);
  PRINT_BLUE
  time_t t = time(nullptr);
  struct tm tm = *localtime(&t);
  fprintf(stderr, "[%02d:%02d:%02d]", tm.tm_hour, tm.tm_min, tm.tm_sec);
  PRINT_RED_B
  fprintf(stderr, "|ERROR!| %s(%d): ", msgHead, line);
  PRINT_DEFAULT
  va_list args;
  va_start(args, what);
  vfprintf(stderr, what, args);
  va_end(args);
  fflush(stderr);
  pthread_mutex_unlock(&print_mx);
}

void debug_printf(const char *funcName, const char *format, ...) {
  if (!print_mx_inited) initDebugMX();
  pthread_mutex_lock(&print_mx);
  PRINT_BLUE
  time_t t = time(nullptr);
  struct tm tm = *localtime(&t);
  printf("[%02d:%02d:%02d]", tm.tm_hour, tm.tm_min, tm.tm_sec);
  PRINT_YELLOW_B
  printf("%s ", funcName);
  PRINT_CYAN
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
  PRINT_DEFAULT
  pthread_mutex_unlock(&print_mx);
}

std::string getFileContentsStr(const std::string &filename) {
  std::ifstream t(filename.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

std::string exec(const char *cmd) {
  std::array<char, 128> buffer{};
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
	throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
	result += buffer.data();
  }
  return result;
}

#define DIFF_TOR 0.001

float calVecDiff(float *A, float *B, unsigned int N) {
  float diff = 0.0;
  bool turnDiff = false;
  size_t start = 0;
  for (int i = 0; i < N; ++i) {
	diff += powf(A[i] - B[i], 2.0);
	if (A[i] - B[i] > DIFF_TOR && !turnDiff) {
	  flushed_printf("\nfind %f != %f, diff[%d, ", A[i], B[i], i);
	  turnDiff = true;
	  start = i;
	}
	if (A[i] - B[i] <= DIFF_TOR && turnDiff) {
	  flushed_printf("%d), distance:%d\n", i, i - start);
	  turnDiff = false;
	}
  }
  if (turnDiff) flushed_printf("%d), distance:%d\n", N, N - start);
  //diff = sqrt(diff)/N;
  return diff;
}

void printMatrix(float *A, size_t height, size_t width) { // print part of matrix if large
  for (int i = 0; i < height && i < 20; ++i) {
	for (size_t j = 0; j < width && j < 20; ++j) {
	  if (A[i * width + j] == 0)printf("\033[1;32m");
	  else printf("\033[0m");
	  flushed_printf("%.2f ", A[i * width + j]);
	}
	flushed_printf("\n");
  }
  flushed_printf("\n");
  PRINT_DEFAULT
}

double getCurrentTime() {
  return omp_get_wtime();
}

__unused double getTime() {
  struct timeval tv{};
  gettimeofday(&tv, nullptr);
  return (double)((int64_t)tv.tv_sec * 1000000 + tv.tv_usec) / 1000000.;
}

void string2Size_t(std::string &line, std::vector<size_t> &vector) {
  std::istringstream ss(line);
  std::string word;
  size_t number;
  do {
	ss >> word;
	if (word == " ") continue;
	try {
	  number = std::stoi(word);
	  vector.push_back(number);
	} catch (std::invalid_argument &err) {}
  } while (ss);
}

void string2Floats(std::string &line, std::vector<float> &vector) {
  std::istringstream ss(line);
  std::string word;
  float number;
  do {
	ss >> word;
	if (word == " ") continue;
	try {
	  number = std::stof(word);
	  vector.push_back(number);
	} catch (std::invalid_argument &err) {}
  } while (ss);
}
