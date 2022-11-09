#include <stdio.h>
#include <string.h>

#include <smmintrin.h>

#include "./stringmatching.h"

// extern "C" const char* sse_strncmp(const char*, int, const char*, int);
// extern "C" const char* sse_strncmp_aligned(const char*, int, const char*, int);

int sse_strcmpeq(volatile const char* needle, volatile const char* haystack){
    volatile size_t i;
    for(i = 0;; i+=16){
        asm goto (
            "movdqu     (%1, %0), %%xmm1 \n"
            "pcmpistri  %3, (%2,%0), %%xmm1 \n"
            "jc         %p5 \n"
            "jz         %p4 \n"
            :  
            : "r"(i),
              "r"(needle),
              "r"(haystack),
              "i"(_SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY)
            : 
               "cc",
               "memory",
               "xmm0"
            : equal, different
        );
    }
equal:
    return 1;
different:
    return 0;
}
