#ifndef BITMAP_HPP
#define BITMAP_HPP

#define WORD_OFFSET(i) ((i) >> 5)
#define BIT_OFFSET(i) ((i) & 0x1f)
#include <cstring>
#include <cassert>

class Bitmap {
public:
    unsigned int size, bmSize;
    unsigned int * data = nullptr;
    unsigned nowSize;
    Bitmap() : size(0), data(NULL) {}
    Bitmap(unsigned int size) : size(size) {
        bmSize = WORD_OFFSET(size)+1;
        data = new unsigned int [bmSize]();
        nowSize = bmSize;
    }
    Bitmap(const Bitmap & t) {
        size = t.size;
        bmSize = t.bmSize;
        nowSize = t.nowSize;
        if(data == nullptr) data = new unsigned int [bmSize];
        for(unsigned i = 0; i < nowSize; i++)
            data[i] = t.data[i];
    }
    ~Bitmap() {
        delete [] data;
    }

    void copy(const Bitmap & t) {
        nowSize = t.nowSize;
        for(unsigned i = 0; i < nowSize; i++)
            data[i] = t.data[i];
    }
    void resize(v_size size_) {
        size = size_;
        if(data != nullptr) delete [] data;
        data = new unsigned int [WORD_OFFSET(size)+1]();
        bmSize = WORD_OFFSET(size)+1;
        nowSize = bmSize;
    }
    void ClearAll() {
        memset(data, 0, sizeof(unsigned)*bmSize);
    }
    void Clear() {
        memset(data, 0, sizeof(unsigned)*nowSize);
    }
    void fill() {
        for(v_size i = 0; i < nowSize; i++) {
            data[i] = 0xffffffff;
        }
    }
    void setNowSize(v_size s) {
        nowSize = WORD_OFFSET(s)+1;
        assert(nowSize <= bmSize);
    }

    bool operator [] (v_size i) {
        return data[WORD_OFFSET(i)] & (1u<<BIT_OFFSET(i));
    }

    bool GetBit(v_size i) {
        return data[WORD_OFFSET(i)] & (1u<<BIT_OFFSET(i));
    }
    void SetBit(v_size i) {
        data[WORD_OFFSET(i)] |= (1u<<BIT_OFFSET(i));
    }
    void ClearBit(v_size i) {
        data[WORD_OFFSET(i)] ^= (1u<<BIT_OFFSET(i));
    }

    void inter(const Bitmap * t) {
        // auto print = [](unsigned*d, unsigned n) {
        //     for(unsigned int i = 0; i < 1; i++) {
        //         for(int j = 0; j < 32; j++) {
        //             if((1u<<j)&d[i]) printf("1");
        //             else printf("0");
        //         }
        //     }
        // };
        // print(data, nowSize);
        // printf("\n");
        // print(t->data, nowSize);
        // printf("\n");
        for(unsigned int i = 0; i < nowSize; i++) {
            data[i] &= t->data[i];
        }
        // print(data, nowSize);
        // printf("\n");
    }

    v_size count() {
        v_size cnt = 0;
        // #pragma omp parallel for reduction(+:cnt)
        for (size_t i = 0; i < nowSize; i++) {
            unsigned int word = data[i], tmp = 0;
            while(word != 0) {
                if (word & 1) tmp++;
                word = word >> 1;
            }
            cnt += tmp;
        }
        return cnt;
    }
};

#endif
//新增clear bit 