#include <smmintrin.h>
#include <stdio.h>
#include <cstring>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <iostream>
#include <bitset>

void CoutM128(__m128i m)
{
    unsigned char *m_ = (unsigned char *) &m;
    for (size_t n = 0; n < 16; ++n) {
        unsigned char c = m_[n];
        for (size_t b = 0; b < 8; ++b) {
            // This shows an actual bit order.
            // If you don't care any bit-level operations,
            // This should be
            //   >> bool ret = c & (1 << (8 - b - 1))
            bool ret = c & (1 << b);
            std::cout << ret;
        }
        std::cout << "|";
    }
    std::cout << std::endl;
}

void CoutM128char(__m128i m)
{
    char *m_ = (char *) &m;
    for (size_t n = 0; n < 16; ++n) {
        // int cast is required
        std::cout << m_[n] << ", ";
    }
    std::cout << std::endl;
}


const char* impl1(const char* needle, int nl, const char* haystack, int hl){
    const __m128i prefix = _mm_loadu_si128(reinterpret_cast<const __m128i*>(needle));
    const __m128i zeros  = _mm_setzero_si128();

    for(size_t i = 0; i < hl; i++){
        const __m128i next8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(haystack+i));
        const __m128i result = _mm_mpsadbw_epu8(next8, prefix, 0);
        const __m128i cmp    = _mm_cmpeq_epi16(result, zeros);
        int mask = _mm_movemask_epi8(cmp) & 0x5555;
        while(mask != 0){
            const int pos = __builtin_ctz(mask) / 2;

            if(memcmp(haystack + i + pos, needle, nl) == 0){
                return haystack + i + pos;
            }

            // unset first set bit
          
            mask ^= 1<<pos;
        }
    }
    return NULL;
}

uint8_t buffer[1024*500 + 1];

extern volatile const char* res = NULL;

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
				volatile char* res = strstr((char*)buffer, s1);
			}
            break;
        
        case 1:
            puts("Impl 1");
            for (i=0; i < iterations; i++)
				res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			break;

        // case 2:
        //     puts("Impl 2");
        //     for (i=0; i < iterations; i++)
		// 		impl2(word.c_str(), word.length(), (char*)buffer, size);
        //     break;

        case 3:
        {
			puts("verify");
			char* libc = strstr((char*)buffer, word.c_str());
			const char* impl1res = impl1(word.c_str(), word.length(), (char*)buffer, size);
			// const char* impl2res = impl2(word.c_str(), word.length(), (char*)buffer, size);
            // const char* naiveres = naive(word.c_str(), word.length(), (char*)buffer, size);
			
			printf("LibC = %p\n", libc);
			printf("Impl1 = %p %s\n",
				impl1res,
				(libc != impl1res) ? "FAILED!!!" : "ok"
            );

            // printf("Impl2 = %p %s\n",
			// 	impl2res,
			// 	(libc != impl2res) ? "FAILED!!!" : "ok"
			// );
            // printf("naive = %p %s\n",
			// 	naiveres,
			// 	(libc != naiveres) ? "FAILED!!!" : "ok"
			// );
				
			if (libc != impl1res)
				return 1;    
            break; 
        }    

        // case 4:
        //     puts("naive");
        //     for (i=0; i < iterations; i++)
		// 		naive(word.c_str(), word.length(), (char*)buffer, size);
		// 	break;

    }

}