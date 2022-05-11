//
// Created by pfxu on 4/9/19.
//

#include "heteroCompLib.h"
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include <unistd.h>
#include <fcntl.h>
#include <linux/fb.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <semaphore.h>

//#include "util.h"
//#include "asm-opt.h"


#define SIZE             (32 * 1024 * 1024)
#define BLOCKSIZE        2048
#ifndef MAXREPEATS
# define MAXREPEATS      10
#endif
#ifndef LATBENCH_COUNT
# define LATBENCH_COUNT  10000000
#endif
#define ALIGN_PADDING    0x100000
#define CACHE_LINE_SIZE  64


const char *KERNEL = "\n \
#define CMD_HELPER(FUNC, NAME) FUNC ## _ ## NAME\n\
#define CMD(FUNC, NAME) CMD_HELPER(FUNC, NAME)\n\
\n\
#define STRIDE (1 << STRIDE_ORDER)\n\
#define GRANULARITY (1 << GRANULARITY_ORDER) \n\
\n\
int reduce_int(int v)    { return v; }\n\
int reduce_int2(int2 v)  { return v.x+v.y; }\n\
int reduce_int4(int4 v)  { return v.x+v.y+v.z+v.w; }\n\
int reduce_int8(int8 v)  { return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7; }\n\
int reduce_int16(int16 v){ return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7+v.s8+v.s9+v.sA+v.sB+v.sC+v.sD+v.sE+v.sF; }\n\
\n\
__kernel void initialize(__global DATATYPE *data, const uint n, const int v) {\n\
	const uint id = get_global_id(0);\n\
	const uint stride = get_global_size(0);\n\
	uint offset = 0;\n\
	while( id+offset<n ){\n\
		data[id+offset] = (DATATYPE)v;\n\
		offset += stride;\n\
	}\n\
}\n\
\n\
__kernel void kernel1(__global DATATYPE *data, const uint n) {\n\
	// Get our global thread ID\n\
	const uint id = get_global_id(0);\n\
	const uint  low_order_id = id & ( (uint)(STRIDE - 1));\n\
	const uint high_order_id = id & (~(uint)(STRIDE - 1));\n\
	const uint index = (high_order_id << GRANULARITY_ORDER) | low_order_id;\n\
	const int localid = get_local_id(0);\n\
	const int group_size = get_local_size(0);\n\
/*	__local uint lcount;\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	if( localid==0 )\n\
		lcount = 0;*/\n\
	// Make sure we do not go out of bounds\n\
//	DATATYPE tmp = (DATATYPE)0;\n\
//	#pragma unroll\n\
	for(int i=0; i<GRANULARITY; i++){\n\
		DATATYPE v = data[index+i*STRIDE];\n\
		int v_sum = CMD(reduce, DATATYPE)(v);\n\
		if( v_sum>0 )\n\
			data[index+i*STRIDE] = (DATATYPE)0;\n\
		//tmp = tmp+data[index+i*STRIDE];\n\
		//data[index+i*STRIDE] = (DATATYPE)(index+i*STRIDE);\n\
	}\n\
/*	if( count )\n\
		atomic_add(&lcount, count);\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	// Atomic reduce in global memory\n\
	if(localid==0 && lcount){\n\
		atomic_add((__global int*)result, lcount);\n\
	}*/\n\
}\n\
";

void aligned_block_copy(int64_t *__restrict dst_,
                        int64_t *__restrict src,
                        int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    while ((size -= 64) >= 0) {
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
    }
}

void aligned_block_copy_backwards(int64_t *__restrict dst_,
                                  int64_t *__restrict src,
                                  int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    src += size / 8 - 1;
    dst += size / 8 - 1;
    while ((size -= 64) >= 0) {
        t1 = *src--;
        t2 = *src--;
        t3 = *src--;
        t4 = *src--;
        *dst-- = t1;
        *dst-- = t2;
        *dst-- = t3;
        *dst-- = t4;
        t1 = *src--;
        t2 = *src--;
        t3 = *src--;
        t4 = *src--;
        *dst-- = t1;
        *dst-- = t2;
        *dst-- = t3;
        *dst-- = t4;
    }
}

/*
 * Walk memory addresses in the backwards direction, but still
 * copy each individual 32 byte block in the forward direction.
 */
void aligned_block_copy_backwards_bs32(int64_t *__restrict dst_,
                                       int64_t *__restrict src,
                                       int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    src += size / 8 - 8;
    dst += size / 8 - 8;
    while ((size -= 64) >= 0) {
        t1 = src[4];
        t2 = src[5];
        t3 = src[6];
        t4 = src[7];
        dst[4] = t1;
        dst[5] = t2;
        dst[6] = t3;
        dst[7] = t4;
        t1 = src[0];
        t2 = src[1];
        t3 = src[2];
        t4 = src[3];
        dst[0] = t1;
        dst[1] = t2;
        dst[2] = t3;
        dst[3] = t4;
        src -= 8;
        dst -= 8;
    }
}

/*
 * Walk memory addresses in the backwards direction, but still
 * copy each individual 64 byte block in the forward direction.
 */
void aligned_block_copy_backwards_bs64(int64_t *__restrict dst_,
                                       int64_t *__restrict src,
                                       int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    src += size / 8 - 8;
    dst += size / 8 - 8;
    while ((size -= 64) >= 0) {
        t1 = src[0];
        t2 = src[1];
        t3 = src[2];
        t4 = src[3];
        dst[0] = t1;
        dst[1] = t2;
        dst[2] = t3;
        dst[3] = t4;
        t1 = src[4];
        t2 = src[5];
        t3 = src[6];
        t4 = src[7];
        dst[4] = t1;
        dst[5] = t2;
        dst[6] = t3;
        dst[7] = t4;
        src -= 8;
        dst -= 8;
    }
}

void aligned_block_copy_pf32(int64_t *__restrict dst_,
                             int64_t *__restrict src,
                             int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    while ((size -= 64) >= 0) {
        __builtin_prefetch(src + 32, 0, 0);
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
        __builtin_prefetch(src + 32, 0, 0);
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
    }
}

void aligned_block_copy_pf64(int64_t *__restrict dst_,
                             int64_t *__restrict src,
                             int size) {
    volatile int64_t *dst = dst_;
    int64_t t1, t2, t3, t4;
    while ((size -= 64) >= 0) {
        __builtin_prefetch(src + 32, 0, 0);
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
        t1 = *src++;
        t2 = *src++;
        t3 = *src++;
        t4 = *src++;
        *dst++ = t1;
        *dst++ = t2;
        *dst++ = t3;
        *dst++ = t4;
    }
}

void aligned_block_fill(int64_t *__restrict dst_,
                        int64_t *__restrict src,
                        int size) {
    volatile int64_t *dst = dst_;
    int64_t data = *src;
    while ((size -= 64) >= 0) {
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
        *dst++ = data;
    }
}

/*
 * Simulate reshuffled memory write accesses to the destination
 * buffer (a kind of "drunken master style" access pattern).
 *
 * See: https://github.com/ssvb/tinymembench/issues/7
 */
void aligned_block_fill_shuffle16(int64_t *__restrict dst_,
                                  int64_t *__restrict src,
                                  int size) {
    volatile int64_t *dst = dst_;
    int64_t data = *src;
    while ((size -= 64) >= 0) {
        dst[0 + 0] = data;
        dst[1 + 0] = data;
        dst[1 + 2] = data;
        dst[0 + 2] = data;
        dst[1 + 4] = data;
        dst[0 + 4] = data;
        dst[0 + 6] = data;
        dst[1 + 6] = data;
        dst += 8;
    }
}

void aligned_block_fill_shuffle32(int64_t *__restrict dst_,
                                  int64_t *__restrict src,
                                  int size) {
    volatile int64_t *dst = dst_;
    int64_t data = *src;
    while ((size -= 64) >= 0) {
        dst[3 + 0] = data;
        dst[0 + 0] = data;
        dst[2 + 0] = data;
        dst[1 + 0] = data;
        dst[3 + 4] = data;
        dst[0 + 4] = data;
        dst[2 + 4] = data;
        dst[1 + 4] = data;
        dst += 8;
    }
}

void aligned_block_fill_shuffle64(int64_t *__restrict dst_,
                                  int64_t *__restrict src,
                                  int size) {
    volatile int64_t *dst = dst_;
    int64_t data = *src;
    while ((size -= 64) >= 0) {
        dst[5] = data;
        dst[2] = data;
        dst[7] = data;
        dst[6] = data;
        dst[1] = data;
        dst[3] = data;
        dst[0] = data;
        dst[4] = data;
        dst += 8;
    }
}


static char *align_up(char *ptr, int align) {
    return (char *) (((uintptr_t) ptr + align - 1) & ~(uintptr_t) (align - 1));
}

void *alloc_four_nonaliased_buffers(void **buf1_, int size1,
                                    void **buf2_, int size2,
                                    void **buf3_, int size3,
                                    void **buf4_, int size4) {
    char **buf1 = (char **) buf1_, **buf2 = (char **) buf2_;
    char **buf3 = (char **) buf3_, **buf4 = (char **) buf4_;
    int antialias_pattern_mask = (ALIGN_PADDING - 1) & ~(CACHE_LINE_SIZE - 1);
    char *buf, *ptr;

    if (!buf1 || size1 < 0)
        size1 = 0;
    if (!buf2 || size2 < 0)
        size2 = 0;
    if (!buf3 || size3 < 0)
        size3 = 0;
    if (!buf4 || size4 < 0)
        size4 = 0;

    ptr = buf =
            (char *) malloc(size1 + size2 + size3 + size4 + 9 * ALIGN_PADDING);
    memset(buf, 0xCC, size1 + size2 + size3 + size4 + 9 * ALIGN_PADDING);

    ptr = align_up(ptr, ALIGN_PADDING);
    if (buf1) {
        *buf1 = ptr + (0xAAAAAAAA & antialias_pattern_mask);
        ptr = align_up(*buf1 + size1, ALIGN_PADDING);
    }
    if (buf2) {
        *buf2 = ptr + (0x55555555 & antialias_pattern_mask);
        ptr = align_up(*buf2 + size2, ALIGN_PADDING);
    }
    if (buf3) {
        *buf3 = ptr + (0xCCCCCCCC & antialias_pattern_mask);
        ptr = align_up(*buf3 + size3, ALIGN_PADDING);
    }
    if (buf4) {
        *buf4 = ptr + (0x33333333 & antialias_pattern_mask);
    }

    return buf;
}

static double bandwidth_bench_helper(int64_t *dstbuf, int64_t *srcbuf,
                                     int64_t *tmpbuf,
                                     int size, int blocksize,
                                     const char *indent_prefix,
                                     int use_tmpbuf,
                                     void (*f)(int64_t *, int64_t *, int),
                                     const char *description) {
    int i, j, loopcount, innerloopcount, n;
    double t1, t2;
    double speed, maxspeed;
    double s, s0, s1, s2;

    /* do up to MAXREPEATS measurements */
    s = s0 = s1 = s2 = 0;
    maxspeed = 0;
    for (n = 0; n < MAXREPEATS; n++) {
        f(dstbuf, srcbuf, size);
        loopcount = 0;
        innerloopcount = 1;
        t1 = gettime();
        do {
            loopcount += innerloopcount;
            if (use_tmpbuf) {
                for (i = 0; i < innerloopcount; i++) {
                    for (j = 0; j < size; j += blocksize) {
                        f(tmpbuf, srcbuf + j / sizeof(int64_t), blocksize);
                        f(dstbuf + j / sizeof(int64_t), tmpbuf, blocksize);
                    }
                }
            } else {
                for (i = 0; i < innerloopcount; i++) {
                    f(dstbuf, srcbuf, size);
                }
            }
            innerloopcount *= 2;
            t2 = gettime();
        } while (t2 - t1 < 0.5);
        speed = (double) size * loopcount / (t2 - t1) / 1000000.;

        s0 += 1;
        s1 += speed;
        s2 += speed * speed;

        if (speed > maxspeed)
            maxspeed = speed;

        if (s0 > 2) {
            s = sqrt((s0 * s2 - s1 * s1) / (s0 * (s0 - 1)));
            if (s < maxspeed / 1000.)
                break;
        }
    }

    if (maxspeed > 0 && s / maxspeed * 100. >= 0.1) {
        printf("%s%-52s : %8.1f MB/s (%.1f%%)\n", indent_prefix, description,
               maxspeed, s / maxspeed * 100.);
    } else {
        printf("%s%-52s : %8.1f MB/s\n", indent_prefix, description, maxspeed);
    }
    return maxspeed;
}

void memcpy_wrapper(int64_t *dst, int64_t *src, int size) {
    memcpy(dst, src, size);
}

void memset_wrapper(int64_t *dst, int64_t *src, int size) {
    memset(dst, src[0], size);
}

typedef struct {
    const char *description;
    int use_tmpbuf;

    void (*f)(int64_t *, int64_t *, int);
} bench_info;

static bench_info warmup[] =
        {
                {"C copy backwards",                  0, aligned_block_copy_backwards},
                {"C copy backwards (32 byte blocks)", 0, aligned_block_copy_backwards_bs32},
                {"C copy backwards (64 byte blocks)", 0, aligned_block_copy_backwards_bs64},
                {nullptr,                             0, nullptr}
        };

static bench_info c_benchmarks[] =
        {
                {"C copy backwards",                         0, aligned_block_copy_backwards},
                {"C copy backwards (32 byte blocks)",        0, aligned_block_copy_backwards_bs32},
                {"C copy backwards (64 byte blocks)",        0, aligned_block_copy_backwards_bs64},
                {"C copy",                                   0, aligned_block_copy},
                {"C copy prefetched (32 bytes step)",        0, aligned_block_copy_pf32},
                {"C copy prefetched (64 bytes step)",        0, aligned_block_copy_pf64},
                {"C 2-pass copy",                            1, aligned_block_copy},
                {"C 2-pass copy prefetched (32 bytes step)", 1, aligned_block_copy_pf32},
                {"C 2-pass copy prefetched (64 bytes step)", 1, aligned_block_copy_pf64},
                {"C fill",                                   0, aligned_block_fill},
                {"C fill (shuffle within 16 byte blocks)",   0, aligned_block_fill_shuffle16},
                {"C fill (shuffle within 32 byte blocks)",   0, aligned_block_fill_shuffle32},
                {"C fill (shuffle within 64 byte blocks)",   0, aligned_block_fill_shuffle64},
                {nullptr,                                    0, nullptr}
        };

static bench_info libc_benchmarks[] =
        {
                {"standard memcpy", 0, memcpy_wrapper},
                {"standard memset", 0, memset_wrapper},
                {nullptr,           0, nullptr}
        };

void bandwidth_bench(int64_t *dstbuf, int64_t *srcbuf, int64_t *tmpbuf,
                     int size, int blocksize, const char *indent_prefix,
                     bench_info *bi) {
    while (bi->f) {
        bandwidth_bench_helper(dstbuf, srcbuf, tmpbuf, size, blocksize,
                               indent_prefix, bi->use_tmpbuf,
                               bi->f,
                               bi->description);
        bi++;
    }
}

static void __attribute__((noinline)) random_read_test(char *zerobuffer,
                                                       int count, int nbits) {
    uint32_t seed = 0;
    uintptr_t addrmask = (1 << nbits) - 1;
    uint32_t v;
    static volatile uint32_t dummy;


#define RANDOM_MEM_ACCESS()                 \
        seed = seed * 1103515245 + 12345;       \
        v = (seed >> 16) & 0xFF;                \
        seed = seed * 1103515245 + 12345;       \
        v |= (seed >> 8) & 0xFF00;              \
        seed = seed * 1103515245 + 12345;       \
        v |= seed & 0x7FFF0000;                 \
        seed |= zerobuffer[v & addrmask];

    while (count >= 16) {
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        count -= 16;
    }
    dummy = seed;
#undef RANDOM_MEM_ACCESS
}

static void __attribute__((noinline)) random_dual_read_test(char *zerobuffer,
                                                            int count, int nbits) {
    uint32_t seed = 0;
    uintptr_t addrmask = (1 << nbits) - 1;
    uint32_t v1, v2;
    static volatile uint32_t dummy;

#define RANDOM_MEM_ACCESS()                 \
        seed = seed * 1103515245 + 12345;       \
        v1 = (seed >> 8) & 0xFF00;              \
        seed = seed * 1103515245 + 12345;       \
        v2 = (seed >> 8) & 0xFF00;              \
        seed = seed * 1103515245 + 12345;       \
        v1 |= seed & 0x7FFF0000;                \
        seed = seed * 1103515245 + 12345;       \
        v2 |= seed & 0x7FFF0000;                \
        seed = seed * 1103515245 + 12345;       \
        v1 |= (seed >> 16) & 0xFF;              \
        v2 |= (seed >> 24);                     \
        v2 &= addrmask;                         \
        v1 ^= v2;                               \
        seed |= zerobuffer[v2];                 \
        seed += zerobuffer[v1 & addrmask];

    while (count >= 16) {
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        RANDOM_MEM_ACCESS();
        count -= 16;
    }
    dummy = seed;
#undef RANDOM_MEM_ACCESS
}

static uint32_t rand32() {
    static int seed = 0;
    uint32_t hi, lo;
    hi = (seed = seed * 1103515245 + 12345) >> 16;
    lo = (seed = seed * 1103515245 + 12345) >> 16;
    return (hi << 16) + lo;
}

int latency_bench(int size, int count, int use_hugepage) {
    double t, t2, t_before, t_after, t_noaccess, t_noaccess2;
    double xs, xs1, xs2;
    double ys, ys1, ys2;
    double min_t, min_t2;
    int nbits, n;
    char *buffer, *buffer_alloc;
    if (posix_memalign((void **) &buffer_alloc, 4 * 1024 * 1024, size) != 0)
        return 0;

    buffer = buffer_alloc;

    if (use_hugepage && madvise(buffer, size, use_hugepage > 0 ?
                                              MADV_HUGEPAGE : MADV_NOHUGEPAGE) != 0) {
        free(buffer_alloc);
        return 0;
    }
    memset(buffer, 0, size);

    for (n = 1; n <= MAXREPEATS; n++) {
        t_before = gettime();
        random_read_test(buffer, count, 1);
        t_after = gettime();
        if (n == 1 || t_after - t_before < t_noaccess)
            t_noaccess = t_after - t_before;

        t_before = gettime();
        random_dual_read_test(buffer, count, 1);
        t_after = gettime();
        if (n == 1 || t_after - t_before < t_noaccess2)
            t_noaccess2 = t_after - t_before;
    }

    printf("\nblock size : single random read / dual random read");
    if (use_hugepage > 0)
        printf(", [MADV_HUGEPAGE]\n");
    else if (use_hugepage < 0)
        printf(", [MADV_NOHUGEPAGE]\n");
    else
        printf("\n");

    for (nbits = 10; (1 << nbits) <= size; nbits++) {
        int testsize = 1 << nbits;
        xs1 = xs2 = ys = ys1 = ys2 = 0;
        for (n = 1; n <= MAXREPEATS; n++) {
            /*
             * Select a random offset in order to mitigate the unpredictability
             * of cache associativity effects when dealing with different
             * physical memory fragmentation (for PIPT caches). We are reporting
             * the "best" measured latency, some offsets may be better than
             * the others.
             */
            int testoffs = (rand32() % (size / testsize)) * testsize;

            t_before = gettime();
            random_read_test(buffer + testoffs, count, nbits);
            t_after = gettime();
            t = t_after - t_before - t_noaccess;
            if (t < 0) t = 0;

            xs1 += t;
            xs2 += t * t;

            if (n == 1 || t < min_t)
                min_t = t;

            t_before = gettime();
            random_dual_read_test(buffer + testoffs, count, nbits);
            t_after = gettime();
            t2 = t_after - t_before - t_noaccess2;
            if (t2 < 0) t2 = 0;

            ys1 += t2;
            ys2 += t2 * t2;

            if (n == 1 || t2 < min_t2)
                min_t2 = t2;

            if (n > 2) {
                xs = sqrt((xs2 * n - xs1 * xs1) / (n * (n - 1)));
                ys = sqrt((ys2 * n - ys1 * ys1) / (n * (n - 1)));
                if (xs < min_t / 1000. && ys < min_t2 / 1000.)
                    break;
            }
        }
        printf("%10d : %6.1f ns          /  %6.1f ns \n", (1 << nbits),
               min_t * 1000000000. / count, min_t2 * 1000000000. / count);
    }
    free(buffer_alloc);
    return 1;
}

//int latency_bench_zero_copy(cl_wrapper *cl, int size, int count) {
int latency_bench_zero_copy(char *ptr, int size, int count) {
    double t, t2, t_before, t_after, t_noaccess, t_noaccess2;
    double xs, xs1, xs2;
    double ys, ys1, ys2;
    double min_t, min_t2;
    int nbits, n;
    char *buffer = ptr;

//    ZeroCopyMem<char> zeroCopyMem = init_zero_copy_region<char>(cl, size);
//    buffer = zeroCopyMem.hostPtr;
    memset(buffer, 0, size);

    for (n = 1; n <= MAXREPEATS; n++) {
        t_before = gettime();
        random_read_test(buffer, count, 1);
        t_after = gettime();
        if (n == 1 || t_after - t_before < t_noaccess)
            t_noaccess = t_after - t_before;

        t_before = gettime();
        random_dual_read_test(buffer, count, 1);
        t_after = gettime();
        if (n == 1 || t_after - t_before < t_noaccess2)
            t_noaccess2 = t_after - t_before;
    }

    printf("\nblock size : single random read / dual random read");
    printf("\n");

    for (nbits = 10; (1 << nbits) <= size; nbits++) {
        int testsize = 1 << nbits;
        xs1 = xs2 = ys = ys1 = ys2 = 0;
        for (n = 1; n <= MAXREPEATS; n++) {
            /*
             * Select a random offset in order to mitigate the unpredictability
             * of cache associativity effects when dealing with different
             * physical memory fragmentation (for PIPT caches). We are reporting
             * the "best" measured latency, some offsets may be better than
             * the others.
             */
            int testoffs = (rand32() % (size / testsize)) * testsize;

            t_before = gettime();
            random_read_test(buffer + testoffs, count, nbits);
            t_after = gettime();
            t = t_after - t_before - t_noaccess;
            if (t < 0) t = 0;

            xs1 += t;
            xs2 += t * t;

            if (n == 1 || t < min_t)
                min_t = t;

            t_before = gettime();
            random_dual_read_test(buffer + testoffs, count, nbits);
            t_after = gettime();
            t2 = t_after - t_before - t_noaccess2;
            if (t2 < 0) t2 = 0;

            ys1 += t2;
            ys2 += t2 * t2;

            if (n == 1 || t2 < min_t2)
                min_t2 = t2;

            if (n > 2) {
                xs = sqrt((xs2 * n - xs1 * xs1) / (n * (n - 1)));
                ys = sqrt((ys2 * n - ys1 * ys1) / (n * (n - 1)));
                if (xs < min_t / 1000. && ys < min_t2 / 1000.)
                    break;
            }
        }
        printf("%10d : %6.1f ns          /  %6.1f ns \n", (1 << nbits),
               min_t * 1000000000. / count, min_t2 * 1000000000. / count);
    }
    //clReleaseMemObject(zeroCopyMem.deviceBuffer);
    return 0;
}

//No memory alignment is faster
void alloc_four_nonaliased_buffers_zero_copy(void *ptr_,
                                             void **buf1_, int size1,
                                             void **buf2_, int size2,
                                             void **buf3_, int size3,
                                             void **buf4_, int size4) {
    char **buf1 = (char **) buf1_, **buf2 = (char **) buf2_;
    char **buf3 = (char **) buf3_, **buf4 = (char **) buf4_;
    char *ptr = (char *) ptr_;

    if (!buf1 || size1 < 0)
        size1 = 0;
    if (!buf2 || size2 < 0)
        size2 = 0;
    if (!buf3 || size3 < 0)
        size3 = 0;
    if (!buf4 || size4 < 0)
        size4 = 0;

    if (buf1) {
        *buf1 = ptr;
        ptr = *buf1 + size1;
    }
    if (buf2) {
        *buf2 = ptr;
        ptr = *buf2 + size2;
    }
    if (buf3) {
        *buf3 = ptr;
        ptr = *buf3 + size3;
    }
    if (buf4) {
        *buf4 = ptr;
    }
}

typedef struct {
    sem_t *semSync;
    cl_command_queue *exeQueue;
    cl_kernel *kernel;
    size_t *glWS;
    size_t *lcWS;
} thread_data_cl;

void *interference_memory_gpu(void *arg) {
    auto *threadData = (thread_data_cl *) arg;

    flushed_printf("GPU interference starts...\n");
    while (true) {
        cl_event ev_wait;
        clEnqueueNDRangeKernel(*threadData->exeQueue, *threadData->kernel, 1,
                               nullptr, threadData->glWS, threadData->lcWS, 0, nullptr, &ev_wait);
        clWaitForEvents(1, &ev_wait);
        if (sem_trywait(threadData->semSync) == 0) break;
    }
    flushed_printf("GPU interference ends...\n");
    pthread_exit(nullptr);
}

typedef struct {
    int64_t *dst;
    int64_t *src;
    int size; // in byte
    int cpu_id;

    void (*f)(int64_t *, int64_t *, int);

    sem_t *syncSem;
    sem_t *mainSyncSem;
} thread_data_cpu;

void *interference_memory(void *arg) {
    auto *threadData = (thread_data_cpu *) arg;
    //if (setCurThreadAffinity(threadData->cpu_id)) {
    //    checkCurThreadAffinity(0);
    //    flushed_printf("pthread launch on core %d failed!\n", threadData->cpu_id);
    //    pthread_exit(nullptr);
    //}

    while (setCurThreadAffinity(threadData->cpu_id) != 0) {//
        //  try to set affinity, a war must win!
    }
    sem_post(threadData->mainSyncSem);

    flushed_printf("Memory interference thread... %d\n", threadData->cpu_id);

    while (true) {
        int64_t *dst_ = threadData->dst;
        int64_t *src_ = threadData->src;
        long size_ = threadData->size;
        threadData->f(dst_, src_, size_);
        if (sem_trywait(threadData->syncSem) == 0) break;
    }

    flushed_printf("Memory interference done...%d\n", threadData->cpu_id);
    pthread_exit(nullptr);
}

int main(int argc, char **argv) {
    if (argc > 1) {
        int affinity = atoi(argv[1]);
        setCurThreadAffinity(affinity);
        checkCurThreadAffinity(0);
    }

    int latbench_size = SIZE * 2, latbench_count = LATBENCH_COUNT;
    int64_t *srcbuf, *dstbuf, *tmpbuf;
    size_t bufsize = SIZE;

    //Init OpenCL & prepare GPU interference
    unsigned int log2_indexes = 26;
    unsigned int log2_grid = 20;
    unsigned int log2_wgroup = 8;
    unsigned int vecsize = 2;
    char s_vecsize[3] = "";
    if (vecsize > 1)
        sprintf(s_vecsize, "%d", vecsize);
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    auto *cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    char options[4096];
    sprintf(options, " -DDATATYPE=int%s -DSTRIDE_ORDER=%d -DGRANULARITY_ORDER=%d",
            s_vecsize, 8, log2_indexes - log2_grid);
    cl_program program = cl->createProgramWithOptions(&KERNEL, 1, options);
    cl_kernel kernel = cl->createKernel("kernel1", program);
    cl_command_queue exeQueue = cl->createProfilingCmdQueue();
    int index_space = pow2(log2_indexes);
    flushed_printf("Applied zero copy memory space: %d MB\n",
                   index_space * vecsize * sizeof(int) / 1024 / 1024);
    auto zeroCopyMem = init_zero_copy_region<int>(cl, index_space * vecsize);
    alloc_four_nonaliased_buffers_zero_copy((void *) zeroCopyMem.hostPtr,
                                            (void **) &srcbuf, bufsize,
                                            (void **) &dstbuf, bufsize,
                                            (void **) &tmpbuf, BLOCKSIZE,
                                            nullptr, 0);

    clSetKernelArg(kernel, 0, sizeof(cl_mem), &zeroCopyMem.deviceBuffer);
    clSetKernelArg(kernel, 1, sizeof(cl_int), &index_space);
    size_t glWS[1] = {index_space / pow2(log2_indexes - log2_grid)};
    size_t lcWS[1] = {pow2(log2_wgroup)};
    int rc;

    flushed_printf("=====No Interference base line, Big core====\n");
    setCurThreadAffinity(7);
    checkCurThreadAffinity(0);
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    flushed_printf("===== Done! ======\n");

    flushed_printf("=====No Interference base line, Little core====\n");
    setCurThreadAffinity(0);
    checkCurThreadAffinity(0);
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    flushed_printf("===== Done! ======\n");

    sem_t syncSem;
    pthread_t gpu_thread;
    thread_data_cl threadData = {
            .semSync  = &syncSem,
            .exeQueue = &exeQueue,
            .kernel   = &kernel,
            .glWS     = glWS,
            .lcWS     = lcWS};
    sem_init(&syncSem, 0, 0);
    if ((rc = pthread_create(&gpu_thread, nullptr, interference_memory_gpu, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);

    flushed_printf("=====GPU interference, Big core====\n");
    setCurThreadAffinity(7);
    checkCurThreadAffinity(0);
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    flushed_printf("===== Done! ======\n");

    flushed_printf("=====GPU interference, Little core====\n");
    setCurThreadAffinity(0);
    checkCurThreadAffinity(0);
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    flushed_printf("===== Done! ======\n");
    sem_post(&syncSem);
    pthread_join(gpu_thread, nullptr);

    pthread_t cpu_threads[8];
    sem_t syncSems[8];
    sem_t mainSyncSems[8];
    thread_data_cpu cpuData[8];

    flushed_printf("=====Multiple CPU interference, Big core====\n");
    setCurThreadAffinity(7);
    checkCurThreadAffinity(0);
    for (int i = 0; i < 7; ++i) {
        sem_init(&syncSems[i], 0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        cpuData[i] = {
                .dst  = dstbuf,
                .src  = srcbuf,
                .size = static_cast<int>(bufsize),
                .cpu_id = i,
                .f = aligned_block_copy_backwards,
                .syncSem = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&cpu_threads[i], nullptr, interference_memory, &cpuData[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 0; i < 7; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    flushed_printf("Interference threads online, benchmark start... \n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    for (int i = 0; i < 7; i++) {
        sem_post(&syncSems[i]);
        pthread_join(cpu_threads[i], nullptr);
    }
    flushed_printf("===== Done! ======\n");

    flushed_printf("=====Multiple CPU + GPU interference, Big core====\n");
    sem_init(&syncSem, 0, 0);
    if ((rc = pthread_create(&gpu_thread, nullptr, interference_memory_gpu, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    setCurThreadAffinity(7);
    checkCurThreadAffinity(0);
    for (int i = 0; i < 7; ++i) {
        sem_init(&syncSems[i], 0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        cpuData[i] = {
                .dst  = dstbuf,
                .src  = srcbuf,
                .size = static_cast<int>(bufsize),
                .cpu_id = i,
                .f = aligned_block_copy_backwards,
                .syncSem = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&cpu_threads[i], nullptr, interference_memory, &cpuData[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 0; i < 7; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    flushed_printf("Interference threads online, benchmark start... \n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    sem_post(&syncSem);
    pthread_join(gpu_thread, nullptr);
    for (int i = 0; i < 7; i++) {
        sem_post(&syncSems[i]);
        pthread_join(cpu_threads[i], nullptr);
    }
    flushed_printf("===== Done! ======\n");


    flushed_printf("=====Multiple CPU interference, Little core====\n");
    setCurThreadAffinity(0);
    checkCurThreadAffinity(0);
    for (int i = 1; i < 8; ++i) {
        sem_init(&syncSems[i], 0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        cpuData[i] = {
                .dst  = dstbuf,
                .src  = srcbuf,
                .size = static_cast<int>(bufsize),
                .cpu_id = i,
                .f = aligned_block_copy_backwards,
                .syncSem = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&cpu_threads[i], nullptr, interference_memory, &cpuData[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 1; i < 8; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    flushed_printf("Interference threads online, benchmark start... \n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    for (int i = 1; i < 8; i++) {
        sem_post(&syncSems[i]);
        pthread_join(cpu_threads[i], nullptr);
    }
    flushed_printf("===== Done! ======\n");

    flushed_printf("=====Multiple CPU + GPU interference, Little core====\n");
    sem_init(&syncSem, 0, 0);
    if ((rc = pthread_create(&gpu_thread, nullptr, interference_memory_gpu, &threadData)))
        fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    setCurThreadAffinity(0);
    checkCurThreadAffinity(0);
    for (int i = 1; i < 8; ++i) {
        sem_init(&syncSems[i], 0, 0);
        sem_init(&mainSyncSems[i], 0, 0);
        cpuData[i] = {
                .dst  = dstbuf,
                .src  = srcbuf,
                .size = static_cast<int>(bufsize),
                .cpu_id = i,
                .f = aligned_block_copy_backwards,
                .syncSem = &syncSems[i],
                .mainSyncSem = &mainSyncSems[i]};
        if ((rc = pthread_create(&cpu_threads[i], nullptr, interference_memory, &cpuData[i])))
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
    }
    for (int i = 1; i < 8; ++i) {
        sem_wait(&mainSyncSems[i]);
    }
    flushed_printf("Interference threads online, benchmark start... \n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", c_benchmarks);
    printf(" ---\n");
    bandwidth_bench(dstbuf, srcbuf, tmpbuf, bufsize, BLOCKSIZE, " ", libc_benchmarks);
    latency_bench_zero_copy((char *) zeroCopyMem.hostPtr, latbench_size, latbench_count);
    sem_post(&syncSem);
    pthread_join(gpu_thread, nullptr);
    for (int i = 1; i < 8; i++) {
        sem_post(&syncSems[i]);
        pthread_join(cpu_threads[i], nullptr);
    }
    flushed_printf("===== Done! ======\n");


    sem_destroy(&syncSem);
    for (int i = 0; i < 8; ++i) {
        sem_destroy(&syncSems[i]);
        sem_destroy(&mainSyncSems[i]);
    }
    clReleaseMemObject(zeroCopyMem.deviceBuffer);
    return 0;
}