#include <smmintrin.h>
#include <stdio.h>
#include <cstring>
#include <cstdint>
#include <string>
#include <unordered_map>

const char* impl1(const char* needle, int nl, const char* haystack, int hl){
    const __m128i first = _mm_set1_epi8(needle[0]);
    const __m128i last = _mm_set1_epi8(needle[nl - 1]);

    const char* haystack_last_start = haystack+nl-1;

    
    for (size_t i = 0; i < hl; i += 16) {

        //loadu is stupid, we can align ourselves once and make it faster?
        //_mm256_loadu_si256
        const __m128i block_first = _mm_loadu_si128(reinterpret_cast<const __m128i *>(haystack+i));
        const __m128i block_last = _mm_loadu_si128(reinterpret_cast<const __m128i *>(haystack_last_start+i));
        const __m128i first_matches = _mm_cmpeq_epi8(block_first, first);
        const __m128i last_matches = _mm_cmpeq_epi8(block_last, last);
        const __m128i ends_matches = _mm_and_si128(first_matches, last_matches);

        // // slightly faster than move mask, requires SSE4
        // const bool isAllZero = _mm_testz_si128(ends_matches, ends_matches);
        int mask = _mm_movemask_epi8(ends_matches);
        while(mask){
            // get first set bit using gcc compiler intrinsics
            const int firstSetBitPos = __builtin_ctz(mask);

            if(memcmp(haystack + i + firstSetBitPos, needle, nl) == 0){
                return haystack + i + firstSetBitPos;
            }

            // unset first set bit
          
            mask ^= 1<<firstSetBitPos;
            // toggling with xor should be faster than clearing with and
            // cause it doesn't take an extra ~ operation to build the mask
            // since we already know the bit is set
            
        }
    }

    return NULL;
}

const char* naive(const char* needle, int nl, const char* haystack, int hl){
    for(int i = 0; i < hl-nl; i++){
        if(haystack[i] == needle[0] && haystack[i+nl-1] == needle[nl-1]){
            if(memcmp(haystack + i, needle, nl) == 0){
                return haystack + i;
            }
        }
    }
    return NULL;
}

const char* impl2(const char* needle, int nl, const char* haystack, int hl){
    const __m128i first = _mm_set1_epi8(needle[0]);
    const __m128i last = _mm_set1_epi8(needle[nl - 1]);

    const char* haystack_last_start = haystack+nl-1;
    if (((intptr_t)haystack & 0xF) != 0) {

        int offset = 16 - ((intptr_t)haystack & 0xF);
        if(offset > nl){
            const char* res = impl1(needle, nl, haystack, offset);
            if(res != NULL)
                return res;
        }
        haystack += offset;
    }

    if (((intptr_t)haystack & 0xF) != 0) {
        printf("NOT ALIGNED!!!");
        return NULL;
    }

    for (size_t i = 0; i < hl; i += 16) {

        //loadu is stupid, we can align ourselves once and make it faster?
        //_mm256_loadu_si256
        
        const __m128i block_first = _mm_load_si128(reinterpret_cast<const __m128i *>(haystack+i));
        const __m128i block_last = _mm_loadu_si128(reinterpret_cast<const __m128i *>(haystack_last_start+i));
        const __m128i first_matches = _mm_cmpeq_epi8(block_first, first);
        const __m128i last_matches = _mm_cmpeq_epi8(block_last, last);
        const __m128i ends_matches = _mm_and_si128(first_matches, last_matches);

        // // slightly faster than move mask, requires SSE4
        // const bool isAllZero = _mm_testz_si128(ends_matches, ends_matches);
        int mask = _mm_movemask_epi8(ends_matches);
        while(mask){
            // get first set bit using gcc compiler intrinsics
            const int firstSetBitPos = __builtin_ctz(mask);

            if(memcmp(haystack + i + firstSetBitPos, needle, nl) == 0){
                return haystack + i + firstSetBitPos;
            }

            // unset first set bit
          
            mask ^= 1<<firstSetBitPos;
            // toggling with xor should be faster than clearing with and
            // cause it doesn't take an extra ~ operation to build the mask
            // since we already know the bit is set
            
        }
    }

    return NULL;
}

uint8_t buffer[1024*500 + 1];

extern volatile char* libcres = NULL;
extern volatile const char* impl1res = NULL;
extern volatile const char* impl2res = NULL;


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
        {"libc", 0},
        {"impl1", 1},
        {"impl2", 2},
        {"verify", 3},
        {"naive", 4}
    })[algorithmVersion];

    const char* s1 = word.c_str();

    switch (alg)
    {
        case 0:
			puts("Lib C");
            for (i=0; i < iterations; i++) {
				libcres = strstr((char*)buffer, s1);
			}
            break;
        
        case 1:
            puts("Impl 1");
            for (i=0; i < iterations; i++)
				impl1res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			break;

        case 2:
            puts("Impl 2");
            for (i=0; i < iterations; i++)
				impl2res = impl2(word.c_str(), word.length(), (char*)buffer, size);
            break;

        case 3:
        {
			puts("verify");
			char* libc = strstr((char*)buffer, word.c_str());
			const char* impl1res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			const char* impl2res = impl2(word.c_str(), word.length(), (char*)buffer, size);
            const char* naiveres = naive(word.c_str(), word.length(), (char*)buffer, size);
			
			printf("LibC = %p\n", libc);
			printf("Impl1 = %p %s\n",
				impl1res,
				(libc != impl1res) ? "FAILED!!!" : "ok"
			);
            printf("Impl2 = %p %s\n",
				impl2res,
				(libc != impl2res) ? "FAILED!!!" : "ok"
			);
            printf("naive = %p %s\n",
				naiveres,
				(libc != naiveres) ? "FAILED!!!" : "ok"
			);
				
			if (libc != impl1res)
				return 1;    
            break; 
        }    

        case 4:
            puts("naive");
            for (i=0; i < iterations; i++)
				naive(word.c_str(), word.length(), (char*)buffer, size);
			break;

    }

}