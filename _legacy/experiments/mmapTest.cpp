//
// Created by pfxu on 28/8/2019.
//
#include "heteroCompLib.h"
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#define PAGE_SIZE (4*1024)
#define START_ADDRESS 0x8000000
#define START_ADDRESS2 0x800000000
#define N_ITERATIONS 131072

// one more iteration and mmap fails
//#define N_ITERATIONS 524289

void allocate(char *base_address) {
    unsigned long i;
    for (i = 0; i < N_ITERATIONS; ++i) {
        char *current_addr;
        current_addr = base_address + PAGE_SIZE * i;
        char *ret = (char *) mmap((char *) current_addr, PAGE_SIZE, PROT_EXEC | PROT_READ | PROT_WRITE,
                                  MAP_NORESERVE | MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS, 0, 0);
        if (ret == MAP_FAILED) {
            fprintf(stderr, "Error mmap. errno: %d\n", errno);
            exit(-1);
        }
//        printf("%lu\n", i);
    }

    memset((void *) START_ADDRESS, 0, PAGE_SIZE * N_ITERATIONS);
}

int main() {
    allocate((char *) START_ADDRESS);
    allocate((char *) START_ADDRESS2);
    return 0;
}

