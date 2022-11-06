#include <smmintrin.h>
#include <stdio.h>
#include <cstring>
#include <cstdint>
#include <string>
#include <unordered_map>

const char* impl1(const char* needle, int nl, const char* haystack, int hl){
    const __m128i N = _mm_loadu_si128((__m128i*)needle);

    for (size_t i = 0; i < hl; i += 16) {

        const int mode = _SIDD_UBYTE_OPS
                       | _SIDD_CMP_EQUAL_ORDERED
                       | _SIDD_BIT_MASK;

        const __m128i D   = _mm_loadu_si128((__m128i*)(haystack + i));
        const __m128i res = _mm_cmpestrm(N, nl, D, hl - i, mode);
        uint64_t mask = _mm_cvtsi128_si64(res);

        while (mask != 0) {

            const int firstSetBitPos = __builtin_ctz(mask);

            if (memcmp(haystack + i + firstSetBitPos, needle , nl) == 0) {
                return haystack + i + firstSetBitPos;
            }

            mask ^= 1<<firstSetBitPos;
        }
    }

    return NULL;
}


uint8_t buffer[1024*500 + 1];

extern volatile const char* impl1res = NULL;


int main(int argc, char** argv){
    std::string filename = std::string(argv[1]);
    std::string algorithmVersion = std::string(argv[2]);
    char* iterations_str = argv[3];
    std::string word = std::string(argv[4]);
    int iterations;
    sscanf(iterations_str, "%d", &iterations);

    FILE* f;
	int i;
	int size;

	f = fopen(filename.c_str(), "r");
	if (!f) {
		printf("can't open '%s'\n", filename.c_str());
		return 2;
	}
		
	size = fread(buffer, 1, sizeof(buffer), f);
	buffer[size] = 0;
	fclose(f);

    int alg = (std::unordered_map<std::string, int> {
        {"impl1", 1},
    })[algorithmVersion];

    const char* s1 = word.c_str();

    switch (alg)
    {
        case 1:
            puts("Impl 1");
            for (i=0; i < iterations; i++)
				impl1res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			break;

        case 3:
        {
			puts("verify");
			char* libc = strstr((char*)buffer, word.c_str());
			const char* impl1res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			
			printf("LibC = %p\n", libc);
			printf("Impl1 = %p %s\n",
				impl1res,
				(libc != impl1res) ? "FAILED!!!" : "ok"
			);
				
			if (libc != impl1res)
				return 1;    
            break; 
        }    
    }

}