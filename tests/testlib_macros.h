#ifndef EDCL_TESTS_TESTLIB_MACROS_H_
#define EDCL_TESTS_TESTLIB_MACROS_H_
#include "utils.h"

#define TIME_CODE_LOOP(TIMER, LOOP_NUM, CODE_BLOCK) \
{ \
    TIMER=getCurrentTime(); \
    for (int loopCnt = 0; loopCnt < LOOP_NUM; ++loopCnt) \
        {CODE_BLOCK} \
    TIMER=(getCurrentTime()-TIMER)/LOOP_NUM; \
}

#define TIME_CODE_ONCE(TIMER, CODE_BLOCK) \
{ \
    TIMER=getCurrentTime(); \
        {CODE_BLOCK} \
    TIMER=(getCurrentTime()-TIMER); \
}



#endif //EDCL_TESTS_TESTLIB_MACROS_H_
