#ifndef UTILS_H
#define UTILS_H

#include<cstdlib>
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include <omp.h>
#include <random>
#include <semaphore.h>

#define MAX(X, Y) X > Y ? X : Y
#define MIN(X, Y) X < Y ? X : Y


unsigned int pow2(unsigned int v);


double sqr(double v);

void flushed_printf(const char *format, ...);

int compare_doubles(const void *a, const void *b);

float calVecDiff(float *A, float *B, int N);

void printMatrix(float *A, size_t height, size_t width);

std::string getFileContentsStr(std::string filename);

double getCurrentTime();

double gettime();

//Random number generate reference: https://channel9.msdn.com/Events/GoingNative/2013/rand-Considered-Harmful
template<typename T>
//T can be int, float, double, etc.
T randNumGenerator(float low, float high);

template<typename T>
void randArrayGenerator(T low, T high, T *arr, size_t len);

//TODO c++ random lib is so slow(snapdragon 845), rand() not really that random. Good random gen needed
template<typename T>
T randNumGenerator(T low, T high) {
    /*******/
    /*
     * Bad old random number generator, but faster
     * */
    return static_cast<T> ((static_cast<float> (rand()) / static_cast<float> (RAND_MAX)) * 1); //cast output
    /**********/
    //return  3.14;
    //return 1.0;
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

//#pragma omp parallel for num_threads(4) private(i)
    for (i = 0; i < len; ++i) {
        arr[i] = randNumGenerator(low, high);
//        flushed_printf("arr[%d] = %f\n", i, arr[i]);
    }
}


template<typename T>
struct ConcurrentQ {
    sem_t lock;
    int front, rear, size;
    unsigned capacity;
    T *array;
};


template<typename T>
struct ConcurrentQ<T> *ConcurrentQ_create(unsigned capacity) {
    auto *queue = (struct ConcurrentQ<T> *) malloc(sizeof(struct ConcurrentQ<T>));
    queue->capacity = capacity;
    queue->front = queue->size = 0;
    queue->rear = capacity - 1;
    queue->array = (T *) malloc(queue->capacity * sizeof(T));
    sem_init(&queue->lock, 0, 1);
    return queue;
}

template<typename T>
void ConcurrentQ_destroy(struct ConcurrentQ<T> *queue) {
    sem_destroy(&queue->lock);
    free(queue->array);
    free(queue);
}

// Queue is full when size becomes equal to the capacity
template<typename T>
bool ConcurrentQ_isFull(struct ConcurrentQ<T> *queue) { return (queue->size == queue->capacity); }

// Queue is empty when size is 0
template<typename T>
bool ConcurrentQ_isEmpty(struct ConcurrentQ<T> *queue) { return (queue->size == 0); }

// Function to add an item to the queue.
// It changes rear and size
template<typename T>
void ConcurrentQ_enqueue(struct ConcurrentQ<T> *queue, T item) {
    sem_wait(&queue->lock);
    if (ConcurrentQ_isFull(queue)) {
        sem_post(&queue->lock);
        return;
    }
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->array[queue->rear] = item;
    queue->size = queue->size + 1;
    sem_post(&queue->lock);
}

// Function to remove an item from queue.
// It changes front and size
template<typename T>
bool ConcurrentQ_dequeueF(struct ConcurrentQ<T> *queue, T *item) {
    sem_wait(&queue->lock);
    if (ConcurrentQ_isEmpty(queue)) {
        sem_post(&queue->lock);
        return false;
    }
    *item = queue->array[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size = queue->size - 1;
    sem_post(&queue->lock);
    return true;
}

template<typename T>
int ConcurrentQ_dequeueN(struct ConcurrentQ<T> *queue, T *item, int N) {
    sem_wait(&queue->lock);
    if (ConcurrentQ_isEmpty(queue)) {
        sem_post(&queue->lock);
        return -1;
    }
    int n = std::min(N, queue->size);
    for (int i = 0; i < n; ++i) {
        item[i] = queue->array[queue->front];
        queue->front = (queue->front + 1) % queue->capacity;
        queue->size = queue->size - 1;
    }
    sem_post(&queue->lock);
    return n;
}


template<typename T>
bool ConcurrentQ_dequeueR(struct ConcurrentQ<T> *queue, T *item) {
    sem_wait(&queue->lock);
    if (ConcurrentQ_isEmpty(queue)) {
        sem_post(&queue->lock);
        return false;
    }
    *item = queue->array[queue->rear];
    queue->rear = (queue->capacity + queue->rear - 1) % queue->capacity;
    queue->size = queue->size - 1;
    sem_post(&queue->lock);
    return true;
}

template<typename T>
bool ConcurrentQ_front(struct ConcurrentQ<T> *queue, T *item) {
    if (ConcurrentQ_isEmpty(queue)) {
        return false;
    }
    item = &queue->array[queue->front];
    return true;
}

template<typename T>
T ConcurrentQ_rear(struct ConcurrentQ<T> *queue, T *item) {
    if (ConcurrentQ_isEmpty(queue)) {
        return false;
    }
    item = &queue->array[queue->rear];
    return true;
}

template<typename T>
T *ConcurrentQ_QtoArray(struct ConcurrentQ<T> *queue) {
    T *arr = malloc(queue->size * sizeof(T));
    for (int i = 0; i < queue->size; ++i) {
        arr[i] = queue->array[(queue->front + i) % queue->capacity];
    }
    return arr;
}

template<typename T>
void ConcurrentQ_print(struct ConcurrentQ<T> *queue) {
    std::cout << "queue: ";
    for (int i = 0; i < queue->size; ++i) {
        std::cout << queue->array[(queue->front + i) % queue->capacity] << " ";
    }
    std::cout << "(size " << queue->size;
    std::cout << ", front " << queue->front;
    std::cout << ", rear " << queue->rear;
    std::cout << ")\n";
}


#endif
