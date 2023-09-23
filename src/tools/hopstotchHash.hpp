#ifndef HOPSTOTCH
#define HOPSTOTCH

#include "../tools/type.hpp"
#ifndef _WIN32
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

#include <immintrin.h>

#include <utility>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#define EMPTY 0xffffffff

constexpr v_size L2_CACHE_LINE_SIZE = 64;
constexpr v_size H = L2_CACHE_LINE_SIZE / sizeof(v_size);

#ifdef __AVX512__
#define setFunction _mm512_set1_epi32
#define cmpFunction _mm512_cmpeq_epi32_mask
#define loadFunction _mm512_loadu_si512
#endif

#ifdef __AVX512__
const __m256i eightEMPTY = _mm256_set1_epi32(EMPTY);
const __m512i sixteenEMPTY = _mm512_set1_epi32(EMPTY);
#endif

class hopstotchHash {
public:
    v_size * v;
    v_size n; //邻接点数量
    v_size roundN, preZero;
    v_size hashKey = 1e9 + 7;
    
public:
#ifndef __SSE__
#define __SSE__
#endif

#ifdef 	__SSE__
#ifndef __AVX512__
    bool containSIMD(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        auto address = reinterpret_cast<const __m256i *>(v + p);
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(address);
        // __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        // if(cmpRes) return true;
        __m256i cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        auto msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif
        eightNextV = _mm256_loadu_si256(address + 1);
        // cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        // if(cmpRes) return true;
        cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif 
    
        return false;
    }
#else
    bool containSIMD(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p));
        __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    
        
        eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p + 8));
        cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    

        return false;
    }
#endif
#endif

    ~hopstotchHash() {
        delete [] v;
    }

    void build(v_size * adjList, v_size n_) {
        n = n_;
        // preZero = __builtin_clz(2 * n - 1);
        // if(n == 0) return;

        if(n > 0) {
            preZero = 0;
            while((1u<<preZero) <= 2*n-1) preZero++;
            preZero = 32 - preZero;

            roundN = 1u << (32 - preZero);
            if(roundN < H) roundN = H;
            v = new v_size[roundN + H];
        }
        else {
            preZero = 32;
            roundN = 0;
            v = new v_size[H];
        }
        
// printf("zero:%u roundN%u, n%u\n", preZero, roundN, n);fflush(stdout);
        memset(v, EMPTY, sizeof(v_size) * (roundN + H));
        for(v_size i = 0; i < n; i++) {
            v_size u = adjList[i];
            if(!insert(u)) {
                bool f = reBuild();
                if(!f) {
                    printf("error build hash table\n");
                    fflush(stdout);
                    exit(-1);
                }
                else insert(u);
            }
        }
        // memcpy(v + roundN, v, sizeof(v_size)*(H-1));
// printf("zero:%u roundN%u\n", preZero, roundN);
        for(v_size i = 0; i < H - 1; i++) {
            v[i+roundN] = v[i];
        }

        for(v_size i = 0; i < n; i++) {
            v_size u = adjList[i];
            if(!contain(u)) {
                printf("error build hash table 2\n");
                fflush(stdout);
                exit(-1);
            }
        }
        // for(v_size i = 0; i < H - 1; i++) {
        //     if(v[i] != v[i + roundN]) {
        //         printf("error st and ed of hashTable\n");fflush(stdout);
        //     }
        // }
    }

    bool insert(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;
        v_size i = 0;
        for(; i < roundN; i++) {
            if(v[p] == EMPTY) break;
            p = (p + 1) % roundN;
            // if(v[p] == u) return true;
        }

        while((p - hashV + roundN) % roundN >= H) {
            v_size t = (p - H + 1 + roundN) % roundN;
            bool f = false;

            for(v_size i = t; i != p; i = (i + 1) % roundN) {
                if((p - hash(v[i]) + roundN) % roundN < H) {
                    v[p] = v[i];
                    v[i] = EMPTY;
                    p = i;
                    f = true;
                    break;
                }
            }

            if(!f) return false;
        }

        v[p] = u;
        return true;
    }

    bool containNormal(v_size u) {
        v_size hashV = hash(u);

        for(v_size i = hashV; i < hashV + H; i++) {
            // v_size j = i % roundN;
            if(v[i] == EMPTY) return false;
            if(v[i] == u) return true;
        }

        return false;
    }

    bool contain(v_size u) {

#ifdef __AVX512__
        return containSIMD512(u);
#else
        return containSIMD(u);
        // return containNormal(u);
#endif
    }
#ifdef __AVX512__
    bool containSIMD512(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m512i sixteenU = _mm512_set1_epi32(u);
        
        // while(true) {
        __m512i sixteenNextV = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(v + p));
        __mmask16 cmpRes = _mm512_cmpeq_epi32_mask(sixteenU, sixteenNextV);
        if(cmpRes) return true;    
        // if(_mm512_cmpeq_epi32_mask(sixteenNextV, sixteenEMPTY)) return false;

        //     p += parallelism;
        // }
        
        return false;
    }
#endif
    bool reBuild() {
        v_size * tmpV = v;
        v = new v_size[roundN];
        memset(v, 0xff, sizeof(v_size) * roundN);

        v_size rebuildTimes = 0;
        while(rebuildTimes < 1000) {
            bool f = true;
            for(v_size i = rebuildTimes; i < roundN + rebuildTimes; i++) {
                v_size j = i % roundN;
                if(tmpV[j] != EMPTY) {
                    if(!insert(tmpV[j])) {
                        f = false;
                        break;
                    }
                }
            }

            if(f) break;
            else {
                rebuildTimes++;
                memset(v, 0xff, sizeof(v_size) * roundN);
            }
        }
        
        delete [] tmpV;

        if(rebuildTimes == 1000) return false;

        return true;
    }
//h(k) = (A∗k mod 2^w) >> (w − r), w = 32, A常数，r=log(roundN)
    v_size hash(v_size t) {
        v_size tmp = hashKey * t;
        v_size ans = tmp >> preZero;
        // assert(ans < roundN);
        return ans;
    }
};



class CuckooHash
{
    using int32 = int32_t;
    const int unfilled = -1;

private:
	/* data */
	int32 capacity;
	int32 mask;
	int32 size;
	int32 buff_size = sizeof(int32);
	int32 *hashtable = nullptr;

	void rehash(int32 **_table) {
		int32 oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		int32 *newhash = new int32[capacity];
		memset((newhash), unfilled, sizeof(int32) * capacity);
		for (int32 i = 0; i < oldcapacity; ++i){
			if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
		}
		std::swap((*_table), newhash);
		delete[] newhash;
	}
	void insert(const int32 &_u, int32 **_table) {
		
		int32 hs = hash1(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		bool use_hash1 = true;
		int32 u = _u;
		for (int32 i = 0; i < mask; ++i) {
			int32 replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			int32 j = 0;
			for (; j < buff_size; ++j) {
				if ((*_table)[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = std::move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++) {
					(*_table)[hs * buff_size + j - 1] =
						std::move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = std::move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = std::move(replaced);
			if (u == unfilled) return;
		}
		rehash(_table);
		insert(u, _table);
	}

	int32 hash1(const int32 x) { return x & mask;}
	int32 hash2(const int32 x) { return ~x & mask;}

public:
	CuckooHash(/* args */) {
		capacity = 0;
		hashtable = nullptr;
		mask = 0;
		size = 0;
	}
	~CuckooHash() {
		if (hashtable != nullptr) {
			delete[] hashtable;
			hashtable = nullptr;
		}
	}

	void reserve(int32 _size) {
		if (capacity >= _size) return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (_size >= mask * buff_size) mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		if (hashtable != nullptr) {
			delete[] hashtable;
			hashtable = nullptr;
		}
		hashtable = new int32[capacity];
		memset(hashtable, unfilled, sizeof(int32) * capacity);
	}

	void insert(const int32 &_u) {
		if (find(_u)) return;
		insert(_u, &hashtable);
		size++;
	}

	bool find(const int32 &_u) {

		int32 hs1 = hash1(_u);

		int32 hs2 = hash2(_u);

// assert(buff_size == 4 && sizeof (int32) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
// if(buff_size*hs1 >= capacity) {
// 	printf("hs1 %d, cap %d\n", hs1, capacity);fflush(stdout);
// }
// assert(buff_size*hs1 < capacity);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
    bool contain(unsigned u) {
        return find(u);
    }
	int32 getcapacity() {return capacity;}
	int32 getmask() {return mask;}
	int32 *gethashtable() {return hashtable;}
};


#endif