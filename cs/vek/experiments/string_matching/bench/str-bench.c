#include <stdio.h>
#include <string.h>

#include "../src/stringmatching.h"

int main(){
    char haystack[] = "The quick brown fox jumps over a lazy dog";
    // char needle1[] = "jumps";
    char needle2[] = "fasdf";
    // const char* resNaive1 = naive(needle1, sizeof needle1 - 1, haystack, sizeof haystack - 1);
    const char* resNaive2 = naive(needle2, sizeof needle2 - 1, haystack, sizeof haystack - 1);
    printf("%p\n", resNaive2);
}