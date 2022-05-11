#pragma ide diagnostic ignored "bugprone-reserved-identifier"
#pragma ide diagnostic ignored "cert-dcl58-cpp"
#pragma ide diagnostic ignored "hicpp-signed-bitwise"
#pragma ide diagnostic ignored "cert-msc50-cpp"
#pragma ide diagnostic ignored "UnusedLocalVariable"
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
#pragma ide diagnostic ignored "OCUnusedMacroInspection"
#ifndef EDCL_UTILS_H
#define EDCL_UTILS_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <cmath>
#include <random>
#include <ctime>
#include <array>

#define MAX(X, Y) (X) > (Y) ? (X) : (Y)
#define MIN(X, Y) (X) < (Y) ? (X) : (Y)

#define PRINT_RED       fprintf(stderr,"\033[0;31m");printf("\033[0;31m");
#define PRINT_RED_B     fprintf(stderr,"\033[1;31m");printf("\033[1;31m");
#define PRINT_GREEN     fprintf(stderr,"\033[0;32m");printf("\033[0;32m");
#define PRINT_GREEN_B   fprintf(stderr,"\033[1;32m");printf("\033[1;32m");
#define PRINT_YELLOW    fprintf(stderr,"\033[0;33m");printf("\033[0;33m");
#define PRINT_YELLOW_B  fprintf(stderr,"\033[1;33m");printf("\033[1;33m");
#define PRINT_BLUE      fprintf(stderr,"\033[0;34m");printf("\033[0;34m");
#define PRINT_BLUE_B    fprintf(stderr,"\033[1;34m");printf("\033[1;34m");
#define PRINT_MAGENTA   fprintf(stderr,"\033[0;35m");printf("\033[0;35m");
#define PRINT_MAGENTA_B fprintf(stderr,"\033[1;35m");printf("\033[1;35m");
#define PRINT_CYAN      fprintf(stderr,"\033[0;36m");printf("\033[0;36m");
#define PRINT_CYAN_B    fprintf(stderr,"\033[1;36m");printf("\033[1;36m");
#define PRINT_DEFAULT   fprintf(stderr,"\033[0m");   printf("\033[0m");

#define CHOOSE(condition, if_true, if_false) condition ? if_true : if_false

#define TEST_CASE(num, what) case num:{what; break;}

#define UNIQUE_PTR_SIZE_T(var, size) std::unique_ptr<size_t []> var = std::make_unique<size_t []>(size); memset(&var[0],0,size*sizeof(size_t))

void string2Floats(std::string &line, std::vector<float> &vector);

void string2Size_t(std::string &line, std::vector<size_t> &vector);

unsigned int pow2(unsigned int v);

__unused double sqr(double v);

void flushed_printf(const char *msg, ...);

void debug_printf(const char *msgHead, const char *msg, ...);

void err_printf(const char *function, int line, const char *what, ...);

__unused int compare_doubles(const void *a, const void *b);

float calVecDiff(float *A, float *B, unsigned int N);

void printMatrix(float *A, size_t height, size_t width);

__unused std::string getFileContentsStr(const std::string &filename);

double getCurrentTime();

__unused double getTime();

//Random number generate reference: https://channel9.msdn.com/Events/GoingNative/2013/rand-Considered-Harmful
template<typename T>
//T can be int, float, double, etc.
T randNumGenerator(float low, float high);

template<typename T>
void randArrayGenerator(T low, T high, T *arr, size_t len);

template<typename T>
T randNumGenerator(T low, T high) {
	/*******/
	/*
	 * Bad old random number generator, but faster
	 * */
	return static_cast<T> (static_cast<float> (rand()) / static_cast<float> (RAND_MAX) * 1); //cast output
	/**********/
}

template<typename T>
T randNumberGeneratorCXX11(T low, T high) {
	/*
	 * Super slow with c++11 method, don't know why
	 * */
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_real_distribution<> distribution(static_cast<int>(low), static_cast<int>(high)); // define the range
	return distribution(eng);
}

template<typename T>
void randArrayGenerator(T low, T high, T *arr, size_t len) {
	int i;
//#pragma omp parallel for num_threads(4) private(i) // enable this phone may die on big array
	for (i = 0; i < len; ++i) {
//	arr[i] = randNumberGeneratorCXX11(low, high); // super slow
		arr[i] = randNumGenerator(low, high);
	}
}

#define SET_BIT(BF, N) BF |= ((u_int16_t) 0x0001 << N)
#define CLR_BIT(BF, N) BF &= ~((u_int16_t) 0x0001 << N)
#define IS_BIT_SET(BF, N) ((BF >> N) & 0x1)

inline void printBits(u_int16_t bits, uint numBits) {
	int i;
	for (i = 0; i < numBits; ++i) {
		flushed_printf(IS_BIT_SET(bits, i) ? "+" : "-");
	}
	flushed_printf("\n");
}

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace std {
template<class T>
struct _Unique_if {
	typedef unique_ptr<T> _Single_object;
};

template<class T>
struct _Unique_if<T[]> {
	typedef unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N>
struct _Unique_if<T[N]> {
	typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args &&... args) {
	return unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
	typedef typename remove_extent<T>::type U;
	return unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args &&...) = delete;
}

template<typename T>
bool isInVector(std::vector<T> &vec, T &t) {
	for (auto &element : vec) {
		if (t == element) return true;
	}
	return false;
}

std::string exec(const char *cmd);

#endif //EDCL_UTILS_H
