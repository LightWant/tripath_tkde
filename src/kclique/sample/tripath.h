#ifndef TRIPATH_H
#define TRIPATH_H

#include "../../graph/graph.hpp"
#include "../../tools/type.hpp"
#include "../../tools/hopstotchHash.hpp"
#include "../../tools/bitmap.hpp"
#include "../../tools/linearSet.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

struct tripath {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    double sumW;
    e_size t = 0;
    e_size sampleTotalTimes = 0;
    double ** dp;
    double * memoryPool = nullptr;

    v_size * pEdge = nullptr;
    v_size * pIdx = nullptr;
    v_size * pColor = nullptr;
    v_size * sortByColor = nullptr;

    //edge graph
    v_size * pEEdge = nullptr;
    v_size * pEIdx = nullptr;

    v_size vCnt, eCnt;

    v_size * clique = nullptr;
    // std::default_random_engine e(rd());
    e_size N = 5000000;

// std::vector<v_size> idClique;
// std::vector<v_size> hittedTimesPerNode;
// std::vector<v_size> cntPathPerNode;

    void init(v_size sz_, std::vector<v_size> & nodes, e_size N_=5000000) {
// idClique.resize(k);
// hittedTimesPerNode.resize(sz_);
// cntPathPerNode.resize(sz_);
        sz = sz_;
        N = N_;
        experiments = new double[sz];

        // auto cmp = [&](v_size a, v_size b) {
        //     return g->color[a] > g->color[b];
        // };
        
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            
            memset(pColor, 0, sizeof(v_size)*(g->cc+1));
            v_size deg = g->pIdx[u+1] - g->pIdx2[u];
            for(v_size i = 0; i < deg; i++) {
                pColor[g->color[g->pEdge[g->pIdx2[u] + i]] + 1]++;
            }

            for(v_size i = 1; i < g->cc; i++) {
                pColor[i] += pColor[i - 1];
            }

            for(v_size i = 0; i < deg; i++) {
                v_size v = g->pEdge[g->pIdx2[u] + i];
                sortByColor[ pColor[g->color[v]]++ ] = v;
            }
// printf("%u\n", i);fflush(stdout);
            memcpy(g->pEdge + g->pIdx2[u], sortByColor, sizeof(v_size)*deg);
// printf("%u\n", i);fflush(stdout);
            double sumD = computeDP(u);

            experiments[i] = sumD;
            sumW += sumD;
        }

        delete [] sortByColor;
    }

    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        sumW = 0.0;
        clique = new v_size[k];

        dp = new double*[g->degeneracy * g->degeneracy];
        memoryPool = new double[g->degeneracy * g->degeneracy * (k+1)]();
        v_size p = 0;
        for(v_size i = 0; i < g->degeneracy * g->degeneracy; i++) {
            dp[i] = memoryPool + p;
            p += k + 1;
        }
        for(v_size i = 0; i < g->degeneracy * g->degeneracy; i++) {
            dp[i][0] = 0;
            dp[i][1] = 1;
        }

        pEdge = new v_size[g->degeneracy*g->degeneracy];
        pIdx = new v_size[g->degeneracy + 1];
        sortByColor = new v_size[g->degeneracy + 1];
        pColor = new v_size[g->cc + 1];
        pEEdge = new v_size[g->degeneracy*g->degeneracy*g->degeneracy];
        pEIdx = new v_size[g->degeneracy*g->degeneracy + 1];
    }

    ~tripath() {
        if(experiments != nullptr) delete [] experiments;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        // if(sortByColor != nullptr) delete [] sortByColor;
        if(clique != nullptr) delete [] clique;
        if(pEEdge != nullptr) delete [] pEEdge;
        if(pEIdx != nullptr) delete [] pEIdx;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
// bool xxx = true;
    double computeDP(v_size u) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        pIdx[0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = sortByColor[i];
            pIdx[i + 1] = pIdx[i]; 

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = sortByColor[j];
// assert(g->color[v] > g->color[w]);
                if(g->color[v] == g->color[w]) continue;
                if(connect(v, w)) {
// assert(g->color[v] > g->color[w]);
                    pEdge[pIdx[i + 1]++] = j;
                }
            }
        }

        pEIdx[0] = 0;
        v_size m = pIdx[outDegree];
        for(v_size u = 0; u < outDegree; u++) {
            for(v_size l = pIdx[u]; l < pIdx[u + 1]; l++) {
                v_size v = pEdge[l];
                pEIdx[l + 1] = pEIdx[l];
                //N(u)\cap N(v)
// if(xxx && l == 0) {
//     printf("%u-%u\n", u, v);
//     for(v_size l = pIdx[u]; l < pIdx[u+1]; l++) {
//         v_size v = pEdge[l];
//         printf("%u ", v);
//     }
//     printf("\n");
//     for(v_size l = pIdx[v]; l < pIdx[v+1]; l++) {
//         v_size v = pEdge[l];
//         printf("%u ", v);
//     }
//     printf("\n");
// }
                v_size i = l + 1, j = pIdx[v];
                while(true) {
                    while(i < pIdx[u + 1] && pEdge[i] < pEdge[j]) i++;
                    if(i == pIdx[u + 1]) break;
                    while(j < pIdx[v + 1] && pEdge[j] < pEdge[i]) j++;
                    if(j == pIdx[v + 1]) break;

                    if(pEdge[i] == pEdge[j]) {
                        pEEdge[pEIdx[l + 1]++] = i;
// if(xxx && l == 0) {
//     printf("%u\n", pEdge[i]);
// }
                        ++i; ++j;
                    }
                }
            }
        }

        for(v_size i = 0; i < m; i++) dp[i][2] = pEIdx[i + 1] - pEIdx[i];
        for(v_size j = 3; j < k; j++) {
            for(v_size i = 0; i < m; i++) {
                dp[i][j] = 0.0;
                for(v_size l = pEIdx[i]; l < pEIdx[i + 1]; l++) {
                    dp[i][j] += dp[pEEdge[l]][j - 1];
                }
            }
        }

        double sumD = 0;
        for(v_size i = 0; i < m; i++) {
            sumD += dp[i][k - 1];
        }
        return sumD;
    }
int nz = 0;
    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d, std::default_random_engine & e) {
        v_size preId = -1;

        double sumD = experiments[id];
// if(sumD != 0) printf("%u ", id);
        double x = d(e);
        
        double sumTmp = 0.0;
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        v_size m = pIdx[outDegree];
        for(v_size u = 0; u < outDegree; u++) {
            for(v_size i = pIdx[u]; i < pIdx[u + 1]; i++) { 
                sumTmp += dp[i][k - 1];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[0] = sortByColor[ u ];
                    clique[1] = sortByColor[ pEdge[i] ];
                    preId = i;
// idClique[0] = u;
// idClique[1] = pEdge[i];
                    break;
                }
            }

            if(preId != v_size(-1)) break;
        }

// if(dp[preId][k - 1] != 0)
// printf("preId %u %u, dp %u\n", preId, pEIdx[preId + 1] - pEIdx[preId], dp[preId][k - 1]);
        for(v_size i = k - 2; i >= 1; i--) {
            sumTmp = sumD = 0.0;
            for(v_size j = pEIdx[preId]; j < pEIdx[preId + 1]; j++) {
                sumD += dp[pEEdge[j]][i];
            }
// bool f = false;
// printf("sumD :%f\n", sumD);
// if(sumD != 0) nz++;
            x = d(e);
            for(v_size j = pEIdx[preId]; j < pEIdx[preId + 1]; j++) {
                sumTmp += dp[pEEdge[j]][i];
// printf("%f-%f ", sumD, sumTmp);
                if(sumTmp + 1e-15 >= x * sumD) {
                    clique[k - i] = sortByColor[ pEdge[ pEEdge[j] ] ];
// idClique[k-i] = pEdge[ pEEdge[j] ];
// printf("ther\n");
// for(v_size l = 0; l < k - i; l++) {
//     if(!connect(clique[k - i], clique[l])) {
//         printf("error %u\n", l);
//     }
// }
                    preId = pEEdge[j];
// f = true;
                    break;
                }
            }
// assert(f);
            
            for(v_size j = 0; j < k - 1 - i; j++) {
                if(!connect(clique[k - i], clique[j])) {
                    return 0;
                }
            }
        }
// for(v_size j = 0; j < k; j++) cntPathPerNode[idClique[j]]++;
    
//         for(v_size i = k-2; i >= 1; i--) {
//             for(v_size j = 0; j < k - 1 - i; j++) {
//                 if(!connect(clique[k - i], clique[j])) {
//                     return 0;
//                 }
//             }
//         }

// for(v_size j = 0; j < k; j++) hittedTimesPerNode[idClique[j]]++;
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes, double expectedN) {
        t = 0;
        sampleTotalTimes = 0;
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        double ans = 0.0;

        // std::discrete_distribution<int> 
        //       distributionForU(experiments, experiments + sz);
        // std::vector<v_size> timesOfU(sz);
        // for(e_size i = 0; i < sampleTimes; i++) {
        //     timesOfU[distributionForU(generator)]++;
        // }

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            e_size expectedSampleTime
                = std::round(sampleTimes * (experiments[i] / sumW) + 1e-15);
            // e_size expectedSampleTime = timesOfU[i];

            if(expectedSampleTime == 0) continue;

            sortByColor = g->pEdge + g->pIdx2[u];
            // sortGraph(u);
            double sumD = computeDP(u);

            v_size tt = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                tt += sampleOneTime(i, u, uiDistribution, generator);
            }
            t += tt;
            ans += 1.0*tt/expectedSampleTime*experiments[i];
            sampleTotalTimes += expectedSampleTime;
        }
        
        if(sampleTotalTimes < 0) {
            std::discrete_distribution<int> 
              udistribution(experiments, experiments + sz);
printf("|not expected %llu \n", sampleTimes - sampleTotalTimes);
              while(sampleTotalTimes < sampleTimes) {
                  int id = udistribution(generator);
                  v_size u = nodes[id];
                  sortByColor = g->pEdge + g->pIdx2[u];
                  // sortGraph(u);
                  computeDP(u);
                  t += sampleOneTime(id, u, uiDistribution, generator);
                  sampleTotalTimes++;
              }
        }
        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
// printf("nz : %d\n", nz);
//         printf("| %.6f %u %u\n", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
        // printf("| %.8f", expectedN / sumW);

// ui maxHitTime = 0;
// for(v_size i = 0; i < sz; i++) {
//     maxHitTime = std::max(maxHitTime, hittedTimesPerNode[i]);
// }
// std::vector<v_size> hitCntTmp(maxHitTime+1);
// for(v_size i = 0; i < sz; i++) {
//     hitCntTmp[hittedTimesPerNode[i]]++;
// }
// printf("maxHitTime:%u\n", maxHitTime);
// for(ui i = 0; i <= maxHitTime; i++) {
//     if(hitCntTmp[i]) printf("hit%u:%u\n", i, hitCntTmp[i]);
// }
// for(v_size i = 0; i < sz; i++) {
//     if(hittedTimesPerNode[i] > 0) {
// // assert(cntPathPerNode[i] > 0);
//         printf(">1hit_%u_%u:%.4f\n", i, nodes[i], 1.0*hittedTimesPerNode[i]/cntPathPerNode[i]);
//     }


// }

        // if(sampleTotalTimes < sampleTimes) {
        //     sampleTotalTimes = sampleTimes;
        // }

        return 1.0 * t / sampleTotalTimes * sumW;
        //return ans;
    }
};

#endif