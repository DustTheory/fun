#include<stdio.h>
#include <smmintrin.h>


extern "C" const char* naive(const char*, int, const char*, int);

// extern "C" const char* sse_strncmp(const char*, int, const char*, int);
// extern "C" const char* sse_strncmp_aligned(const char*, int, const char*, int);

const char* sse_strcmpeq(const char* needle, int nl, const char* haystack, int hl){
    for(size_t i = 0; i < hl; i+=16){
    //     movdqu     xmm1, [rdi + rdx]
	// pcmpistri  xmm1, [rsi + rdx], 0x18	; EQUAL_EACH | NEGATIVE_POLARITY
	// jc         .dif
	// jz         .eql
    }
}

int main(){
    char haystack[] = "The quick brown fox jumps over a lazy dog";
    char needle[] = "fox";
    const char* location = naive(needle, sizeof needle - 1, haystack, sizeof haystack -1);
    if(location == NULL){
        printf("not found\n");
    } else {
        int index = (int)(location - haystack);
        printf("%d\n", index);
    }
    return 0;
}