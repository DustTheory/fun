#include <stdio.h>
#include <string.h>
#include <smmintrin.h>
#include <benchmark/benchmark.h>


extern const char* naive(const char*, int, const char*, int);

// extern "C" const char* sse_strncmp(const char*, int, const char*, int);
// extern "C" const char* sse_strncmp_aligned(const char*, int, const char*, int);

int sse_strcmpeq(const char* needle, char* haystack){
    size_t i;
    for(i = 0;; i+=16){
        asm goto (
            "movdqu     (%1, %0), %%xmm0 \n"
            "pcmpistri  $0x18, (%2,%0), %%xmm0 \n"
            "jc         %p4 \n"
            "jz         %p3 \n"
            :   
            : "r"(i),
              "r"(needle),
              "r"(haystack)
            : 
               "cc",
               "memory"
            : equal, different
        );
    }
equal:
    return 1;
different:
    return 0;
}

inline int native(char* str1, char* str2){
    return strcmp(str1, str2) == 0;
}

static void native_test1(benchmark::State& state) {

  char str1[] = "Carroline";
  char str2[] = "Carroline";
  for (auto _ : state) {
    native(str1, str2);
  }
}

BENCHMARK(native_test1);

static void unaligned_test1(benchmark::State& state) {

  char str1[] = "Carroline";
  char str2[] = "Carroline";
  for (auto _ : state) {
    sse_strcmpeq(str1, str2);
  }
}

BENCHMARK(unaligned_test1);


int main(){
    // char haystack[] = "The quick brown fox jumps over a lazy dog";
    // char needle[] = "fox";
    // const char* location = naive(needle, sizeof needle - 1, haystack, sizeof haystack -1);
    // if(location == NULL){
    //     printf("not found\n");
    // } else {
    //     int index = (int)(location - haystack);
    //     printf("%d\n", index);
    // }
    char str1[] = "Carroline";
    char str2[] = "Carroline";
    int res = sse_strcmpeq(str1, str2);

    printf("%d\n", res);
    return 0;
}