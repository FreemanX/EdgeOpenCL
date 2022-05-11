#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <locale.h>

#define NUMBYTES (1024*1024*1024)

char bytes[NUMBYTES];
int  batch_size = BATCH_SIZE;
int  cache_line_size = CACHE_LINE_SIZE;

double get_time_in_seconds(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + 1.e-6 * (double)tv.tv_usec;
}

void main(void) {
    int i;
    int j;
    int p_start[batch_size];
    int p[batch_size];
    uint64_t c[batch_size];
    int incs[batch_size];
    int incs_total = 0;
    setlocale(LC_NUMERIC, "");
#ifdef CACHE_LINE_FRIENDLY
    int inc = 1 + (cache_line_size * 0);
#endif
#ifdef CACHE_LINE_FRIENDLY_ISH
    int inc = 1 + (cache_line_size * 1);
#endif
#ifdef CACHE_LINE_FRIENDLIER
    int inc = 1 + (cache_line_size * 8);
#endif
#ifdef CACHE_LINE_UNFRIENDLY
    int inc = 1 + (cache_line_size * 8191);
#endif
#ifdef CACHE_LINE_PREFETCH
    char prefetch_text[] = "with    prefetch";
#else
    char prefetch_text[] = "without prefetch";
#endif
    for(i = 0; i < NUMBYTES; i++) {
        bytes[i] = i & 255;
    }
    for(i = 0; i < batch_size; i++) {
        p_start[i] = ((NUMBYTES / batch_size) * i) + (rand() & 8191);
        p[i] = p_start[i];
        c[i] = 0;
        incs[i] = 0;
    }
    double t1 = get_time_in_seconds();
    for(j = 0; j < (500000000 / batch_size); j++) {
#ifdef CACHE_LINE_PREFETCH
        for(i = 0; i < batch_size; i++) {
            __builtin_prefetch(&bytes[p[i]], 1, 3); // https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
        }
#endif
        for(i = 0; i < batch_size; i++) {
            c[i] += bytes[p[i]];
            incs[i] ++;
            incs_total ++;
            p[i] += inc;
            p[i] = p[i] & (NUMBYTES - 1);
        }
    }
    double t2 = get_time_in_seconds();
    //for(i = 0; i < batch_size; i++) {
    //    printf("- p[%2u]: grand total sum of all bytes with detected cache line size %u at p_start %'13u with %'11u incs of size %'7u %s: %'lu\n", i, cache_line_size, p_start[i], incs[i], inc, &prefetch_text[0], c[i]);
    //}
    
    double t3 = get_time_in_seconds();
    for(j = 0; j < 500000000 ; j++) { //TODO: make this equivalant to bached version 
	    c[j] += bytes[p[j]];
	    incs[j] ++;
	    p[j] += inc;
	    p[j] = p[j] & (NUMBYTES - 1);
    }
    double t4 = get_time_in_seconds();

    printf("- %'u incs in %f sec, base line %f sec, %s using batch_size %3u and inc %'7u\n", incs_total, t2 - t1, t4 - t3, &prefetch_text[0], batch_size, inc);
    //printf("- %'u incs in %f seconds or %'13.0f incs per second %s using batch_size %3u and inc %'7u\n", incs_total, t2 - t1, incs_total / (t2 - t1), &prefetch_text[0], batch_size, inc);
}
