#ifndef CCCPATH_H
#define CCCPAHT_H

#include "../../graph/graph.hpp"
#include "../../tools/type.hpp"
#include "../../tools/hopstotchHash.hpp"
#include "../../tools/bitmap.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

// constexpr v_size batchSize = 50;

struct ccpathLocalColorSampler {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    double sumW;
    double ** dp;
    double * memoryPool = nullptr;

    v_size * pEdge = nullptr;
    v_size * pIdx = nullptr;
    v_size * pColor = nullptr;
    v_size * sortByColor = nullptr;
    v_size vCnt, eCnt;

    v_size * clique = nullptr;
    std::default_random_engine e;
    e_size N = 5000000;

    //local color
    v_size sumDeg = 0;
    v_size * color = nullptr;//e
    v_size * maxColor = 0;//v

    void init(v_size sz_, std::vector<v_size> & nodes, e_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];

        //local color
        sumDeg = 0;
        v_size maxDeg = 0;
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u+1] - g->pIdx2[u];
            sumDeg += deg;
            maxDeg = std::max(maxDeg, deg);
        }

        color = new v_size[sumDeg + 1]();
        maxColor = new v_size[sz]();
        v_size * f = new v_size[maxDeg+1]();
        v_size p = 0;
        v_size cc = 0;

        // for(v_size i = 0; i < sz; i++) {
        //     v_size u = nodes[i];
        //     v_size deg = g->pIdx[u+1] - g->pIdx2[u];

        //     for(v_size x = g->pIdx2[u]; x < g->pIdx[u + 1]; x++) {
        //         v_size v = g->pEdge[x];

        //         for(v_size y = g->pIdx2[u]; y < x; y++) {
        //             v_size w = g->pEdge[y];

        //             if(connect(w, v)) {
        //                 f[ color[ p + y - g->pIdx2[u] ] ] = v;
        //             }
        //         }

        //         v_size c = 0;
        //         while(c+1 < deg && f[c] == v) c++;

        //         color[p + x - g->pIdx2[u]] = c;
        //         maxColor[i] = std::max(maxColor[i], c);
        //     }

        //     cc = std::max(cc, maxColor[i]);
        //     p += deg;
        // }

        std::vector<v_size> degree2hop(maxDeg);
        std::vector<v_size> sortByDeg(maxDeg);
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u+1] - g->pIdx2[u];

            std::fill(degree2hop.begin(), degree2hop.begin() + deg, 0);
            for(v_size x = g->pIdx2[u]; x < g->pIdx[u + 1]; x++) {
                v_size v = g->pEdge[x];
                sortByDeg[x - g->pIdx2[u]] = x - g->pIdx2[u];
            
                for(v_size y = x + 1; y < g->pIdx[u + 1]; y++) {
                    v_size w = g->pEdge[y];

                    if(connect(w, v)) {
                        degree2hop[x - g->pIdx2[u]]++;
                        degree2hop[y - g->pIdx2[u]]++;
                    }
                }
            }
            std::sort(sortByDeg.begin(), sortByDeg.begin() + deg, [&](v_size a, v_size b){
                return degree2hop[a] < degree2hop[b];
            });
// printf("there\n");
            color[p + sortByDeg[0]] = 0;
            for(v_size xx = 1; xx < deg; xx++) {
                v_size x = sortByDeg[xx];
                v_size v = g->pEdge[x + g->pIdx2[u]];

                for(v_size yy = 0; yy < xx; yy++) {
                    v_size y = sortByDeg[yy];
                    v_size w = g->pEdge[y + g->pIdx2[u]];
                    if(connect(v, w)) {
                        f[color[p + y]] = v;
                    }
                }

                v_size c = 0;
                while(c + 1 < deg && f[c] == v) c++;
                color[p + x] = c;
                maxColor[i] = std::max(maxColor[i], c);
            }

            cc = std::max(cc, maxColor[i]);
            p += deg;
        }

        printf("|cc%d ", cc);

        assert(p == sumDeg);

        //sort by local color
        for(v_size i = 0, p = 0; i < sz; i++) {
            v_size u = nodes[i];
            
            memset(pColor, 0, sizeof(v_size)*(maxColor[i] + 1));
            v_size deg = g->pIdx[u+1] - g->pIdx2[u];

            for(v_size j = 0; j < deg; j++) {
                pColor[color[p + j] + 1]++;
            }

            for(v_size j = 1; j < maxColor[i] + 1; j++) {
                pColor[j] += pColor[j - 1];
            }

            for(v_size j = 0; j < deg; j++) {
                v_size v = g->pEdge[g->pIdx2[u] + j];
                sortByColor[ pColor[color[p + j]]++ ] = v;
            }

            memcpy(g->pEdge + g->pIdx2[u], sortByColor, sizeof(v_size)*deg);

            computeDP(u);

            double sumD = 0.0;
            for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                sumD += dp[i][k];
            }

            experiments[i] = sumD;
            sumW += sumD;
            p += deg;
        }

        delete [] sortByColor;
        delete [] f;
    }

    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        sumW = 0.0;
        clique = new v_size[k];

        dp = new double*[g->degeneracy];
        memoryPool = new double[g->degeneracy * (k+1)]();
        v_size p = 0;
        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i] = memoryPool + p;
            p += k + 1;
        }

        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i][0] = 0;
            dp[i][1] = 1;
        }

        pEdge = new v_size[g->degeneracy*g->degeneracy];
        pIdx = new v_size[g->degeneracy + 1];
        sortByColor = new v_size[g->degeneracy + 1];
        pColor = new v_size[g->degeneracy + 1];
    }


    ~ccpathLocalColorSampler() {
        if(experiments != nullptr) delete [] experiments;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        // if(sortByColor != nullptr) delete [] sortByColor;
        if(clique != nullptr) delete [] clique;
        if(color != nullptr) delete [] color;
        if(maxColor != nullptr) delete [] maxColor;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
    
    void computeDP(v_size u) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        pIdx[0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = sortByColor[i];
            pIdx[i + 1] = pIdx[i]; 

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = sortByColor[j];

                // if(color[] == color[]) continue;
                if(hashTable[v].contain(w)) {
                    pEdge[pIdx[i + 1]++] = j;
                }
            }
        }

        for(v_size j = 2; j <= k; j++) {
            for(v_size i = 0; i < outDegree; i++) {
                dp[i][j] = 0.0;
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    dp[i][j] += dp[pEdge[l]][j - 1];
                }
            }
        }
    }

    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d) {
        v_size preId = -1;

        double sumD = experiments[id];
        double x = d(e);
        
        double sumTmp = 0.0;
        v_size deg = g->pIdx[u+1] - g->pIdx2[u];
        for(v_size i = 0; i < deg; i++) {
            sumTmp += dp[i][k];
            if(sumTmp + 1e-10 >= x * sumD) {
                clique[0] = sortByColor[ i ];
                preId = i;
                break;
            }
        }

        if(preId == (v_size)4294967295U) {
// printf("preId %lld|%f,", preId, sumD);fflush(stdout);

            return 0;
        }

        for(v_size i = 1; i < k; i++) {
            sumTmp = sumD = 0.0;
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) {
                sumD += dp[pEdge[j]][k - i];
            }
// bool f = false;
            x = d(e);
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) {
                sumTmp += dp[pEdge[j]][k - i];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[i] = sortByColor[ pEdge[j] ];
                    preId = pEdge[j];
// f = true;
                    break;
                }
            }
// assert(f);
            for(v_size j = 0; j < i-1; j++) {
                if(!connect(clique[i], clique[j])) {
                    return 0;
                }
            }
        }
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes, double expectedN) {
        e_size t = 0;
        e_size sampleTotalTimes = 0;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        double ans = 0.0;

        for(v_size i = 0; i < sz; i++) {
// printf("%ld ", i);fflush(stdout);
            v_size u = nodes[i];

            e_size expectedSampleTime
                = std::round(sampleTimes * (experiments[i] / sumW) + 0.000001);

            if(expectedSampleTime == 0) continue;

            sortByColor = g->pEdge + g->pIdx2[u];
            // sortGraph(u);
            computeDP(u);
// printf("%ld ", i);fflush(stdout);
            v_size tt = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                tt += sampleOneTime(i, u, uiDistribution);
            }
// printf("%ld\n", i);fflush(stdout);
            t += tt;
            ans += 1.0*tt/expectedSampleTime*experiments[i];
            sampleTotalTimes += expectedSampleTime;
        }
        
        if(sampleTotalTimes < sampleTimes) {
            std::discrete_distribution<int> 
              udistribution(experiments, experiments + sz);
printf("|not expected %llu ", sampleTimes - sampleTotalTimes);
              while(sampleTotalTimes < sampleTimes) {
                  int id = udistribution(generator);
                  v_size u = nodes[id];
                  sortByColor = g->pEdge + g->pIdx2[u];
                  // sortGraph(u);
                  computeDP(u);
                  t += sampleOneTime(id, u, uiDistribution);
                  sampleTotalTimes++;
              }
             // printf("|small %.6f %u %u", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
             // return 1.0 * t / sampleTotalTimes * sumW;
        }
        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
        printf("| %.6f %u %u", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
        // printf("| %.8f", expectedN / sumW);
        return 1.0 * t / sampleTotalTimes * sumW;
        // return ans;
    }
};

#endif