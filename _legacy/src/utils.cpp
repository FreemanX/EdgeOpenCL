
#include "utils.h"

unsigned int pow2(unsigned int v) {
    return static_cast<unsigned int>(1 << v);
}

double sqr(double v) {
    return v * v;
}

int compare_doubles(const void *a, const void *b) {
    const auto *v1 = (const double *) a;
    const auto *v2 = (const double *) b;
    if (*v1 == *v2) return 0;
    return *v1 > *v2 ? 1 : -1;
}

void flushed_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    fflush(stdout);
}

std::string getFileContentsStr(std::string filename) {
    std::ifstream t(filename.c_str());
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

#define DIFF_TOR 0.001

float calVecDiff(float *A, float *B, int N) {
    float diff = 0.0;
    bool turnDiff = false;
    size_t start = 0;
    for (int i = 0; i < N; ++i) {
        diff += pow(A[i] - B[i], 2.0);
        if (A[i] - B[i] > DIFF_TOR && !turnDiff) {
            flushed_printf("find %f != %f, diff[%d, ", A[i], B[i], i);
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
    for (int i = 0; i < height && i < 15; ++i) {
        for (size_t j = 0; j < width && j < 15; ++j) {
            //std::cout << A[i * width + j] << " ";
            flushed_printf("%.2f ", A[i * width + j]);
        }
        //std::cout << std::endl;
        flushed_printf("\n");
    }
    //std::cout << std::endl;
    flushed_printf("\n");
}

double getCurrentTime() {
    return omp_get_wtime();
}

//double getCurrentTime() {
//    struct timeval tv{};
//    gettimeofday(&tv, nullptr);
//    return (double) ((int64_t) tv.tv_sec * 1000000 + tv.tv_usec) / 1000000.;
//}

double gettime() {
    struct timeval tv{};
    gettimeofday(&tv, nullptr);
    return (double) ((int64_t) tv.tv_sec * 1000000 + tv.tv_usec) / 1000000.;
}


