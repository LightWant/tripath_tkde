#ifndef TRIPATHPARALLEL_H
#define TRIPATHPARALLEL_H

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

struct tripathParallel {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    v_size threads;

    double * experiments;
    double sumW;
    double *** dp;
    double * memoryPool = nullptr;

    v_size ** pEdge = nullptr;
    v_size ** pIdx = nullptr;
    v_size ** pColor = nullptr;
    v_size ** sortByColor = nullptr;

    //edge graph
    v_size ** pEEdge = nullptr;
    v_size ** pEIdx = nullptr;

    v_size vCnt, eCnt;

    v_size ** clique = nullptr;

    e_size N = 5000000;

    v_size chunkSize = 16;

    void init(v_size sz_, std::vector<v_size> & nodes, e_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];

        #pragma omp parallel reduction(+:sumW)
        {
            int threadId = omp_get_thread_num();

            #pragma omp for schedule(dynamic, chunkSize) 
            for(v_size i = 0; i < sz; i++) {
                v_size u = nodes[i];
                memset(pColor[threadId], 0, sizeof(v_size)*(g->cc+1));
                v_size deg = g->pIdx[u+1] - g->pIdx2[u];
                for(v_size i = 0; i < deg; i++) {
                    pColor[threadId][ g->color[g->pEdge[g->pIdx2[u] + i]] + 1 ]++;
                }

                for(v_size i = 1; i < g->cc; i++) {
                    pColor[threadId][i] += pColor[threadId][i - 1];
                }

                for(v_size i = 0; i < deg; i++) {
                    v_size v = g->pEdge[g->pIdx2[u] + i];
                    sortByColor[threadId][ pColor[threadId][g->color[v]]++ ] = v;
                }

                memcpy(g->pEdge + g->pIdx2[u], sortByColor[threadId], sizeof(v_size)*deg);

                double sumD = computeDP(u, threadId);

                experiments[i] = sumD;
                sumW += sumD;
            }

            delete [] sortByColor[threadId];
            delete [] pColor[threadId];
        }
    }

    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_, v_size threads_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        threads = threads_;

        sumW = 0.0;
        clique = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            clique[i] = new v_size[k];
        }

        dp = new double**[threads];
        for(v_size t = 0; t < threads; t++) {
            // double * memoryPool = new double[g->degeneracy * (k+1)]();
            dp[t] = new double*[g->degeneracy * g->degeneracy];
            // v_size p = 0;

            for(v_size i = 0; i < g->degeneracy * g->degeneracy; i++) {
                dp[t][i] = new double[k+1];
                // p += k + 1;
            }
        }

        for(v_size t = 0; t < threads; t++) {
            for(v_size i = 0; i < g->degeneracy * g->degeneracy; i++) {
                dp[t][i][0] = 0;
                dp[t][i][1] = 1;
            }
        }

        pEdge = new v_size*[threads];
        pIdx = new v_size*[threads];
        for(v_size t = 0; t < threads; t++) {
            pEdge[t] = new v_size[g->degeneracy*g->degeneracy];
            pIdx[t] = new v_size[g->degeneracy + 1];
        }
        sortByColor = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            sortByColor[i] = new v_size[g->degeneracy + 1];
        }
        pColor = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            pColor[i] = new v_size[g->cc + 1];
        }

        pEEdge = new v_size * [threads];
        pEIdx = new v_size * [threads];
        for(v_size i = 0; i < threads; i++) {
            pEEdge[i] = new v_size[g->degeneracy*g->degeneracy*g->degeneracy];
            pEIdx[i] = new v_size[g->degeneracy*g->degeneracy + 1];
        }
    }

    ~tripathParallel() {
        if(experiments != nullptr) delete [] experiments;
        for(v_size i = 0; i < threads; i++) {
            for(v_size j = 0; j < g->degeneracy; j++)
                delete [] dp[i][j];
            delete [] dp[i];
        }
        if(dp != nullptr) delete [] dp;

        for(v_size i = 0; i < threads; i++) {
            delete [] pEdge[i];
            delete [] pIdx[i];
        }
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;

        for(v_size i = 0; i < threads; i++) {
            delete [] clique[i];
        }
        if(clique != nullptr) delete [] clique;

        for(v_size i = 0; i < threads; i++) {
            delete [] pEEdge[i];
            delete [] pEIdx[i];
        }
        delete [] pEEdge;
        delete [] pEIdx;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    double computeDP(v_size u, v_size tId) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        pIdx[tId][0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = sortByColor[tId][i];
            pIdx[tId][i + 1] = pIdx[tId][i]; 

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = sortByColor[tId][j];

                if(g->color[v] == g->color[w]) continue;
                if(connect(v, w)) {
                    pEdge[tId][pIdx[tId][i + 1]++] = j;
                }
            }
        }

        pEIdx[tId][0] = 0;
        v_size m = pIdx[tId][outDegree];
        for(v_size u = 0; u < outDegree; u++) {
            for(v_size l = pIdx[tId][u]; l < pIdx[tId][u + 1]; l++) {
                v_size v = pEdge[tId][l];
                pEIdx[tId][l + 1] = pEIdx[tId][l];

                v_size i = l + 1, j = pIdx[tId][v];
                while(true) {
                    while(i < pIdx[tId][u + 1] && pEdge[tId][i] < pEdge[tId][j]) i++;
                    if(i == pIdx[tId][u + 1]) break;
                    while(j < pIdx[tId][v + 1] && pEdge[tId][j] < pEdge[tId][i]) j++;
                    if(j == pIdx[tId][v + 1]) break;

                    if(pEdge[tId][i] == pEdge[tId][j]) {
                        pEEdge[tId][pEIdx[tId][l + 1]++] = i;
                        ++i; ++j;
                    }
                }
            }
        }

        for(v_size i = 0; i < m; i++) dp[tId][i][2] = pEIdx[tId][i + 1] - pEIdx[tId][i];
        for(v_size j = 3; j < k; j++) {
            for(v_size i = 0; i < m; i++) {
                dp[tId][i][j] = 0.0;
                for(v_size l = pEIdx[tId][i]; l < pEIdx[tId][i + 1]; l++) {
                    dp[tId][i][j] += dp[tId][pEEdge[tId][l]][j - 1];
                }
            }
        }

        double sumD = 0;
        for(v_size i = 0; i < m; i++) {
            sumD += dp[tId][i][k - 1];
        }
        return sumD;
    }

    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d, std::default_random_engine & e, v_size tId) {
        v_size preId = -1;

        double sumD = experiments[id];
        double x = d(e);
        
        double sumTmp = 0.0;
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
      //  v_size m = pIdx[tId][outDegree];
        for(v_size u = 0; u < outDegree; u++) {
            for(v_size i = pIdx[tId][u]; i < pIdx[tId][u + 1]; i++) { 
                sumTmp += dp[tId][i][k - 1];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[tId][0] = sortByColor[tId][ u ];
                    clique[tId][1] = sortByColor[tId][ pEdge[tId][i] ];
                    preId = i;
                    break;
                }
            }

            if(preId != v_size(-1)) break;
        }

        for(v_size i = k - 2; i >= 1; i--) {
            sumTmp = sumD = 0.0;
            for(v_size j = pEIdx[tId][preId]; j < pEIdx[tId][preId + 1]; j++) {
                sumD += dp[tId][pEEdge[tId][j]][i];
            }

            x = d(e);
            for(v_size j = pEIdx[tId][preId]; j < pEIdx[tId][preId + 1]; j++) {
                sumTmp += dp[tId][pEEdge[tId][j]][i];

                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[tId][k - i] = sortByColor[tId][ pEdge[tId][ pEEdge[tId][j] ] ];

                    preId = pEEdge[tId][j];
                    break;
                }
            }
            
            for(v_size j = 0; j < k - 1 - i; j++) {
                if(!connect(clique[tId][k - i], clique[tId][j])) {
                    return 0;
                }
            }
        }
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes, double expectedN) {
        e_size t = 0;
        e_size sampleTotalTimes = 0;
        std::random_device rd;
        // std::default_random_engine generator(rd());
        std::default_random_engine e[200];
        for(v_size i = 0; i < threads; i++) {
            e[i].seed(rd());
        }

        e_size * chunk = new e_size[sz]();
        v_size * pChunk = new v_size[threads + 1];
        pChunk[0] = 0;

        e_size sumCompute = 0;
        #pragma omp parallel for schedule(dynamic, 32) reduction(+:sumCompute, sampleTotalTimes )
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            e_size expectedTime = std::round(sampleTimes * (experiments[i] / sumW));
            if(expectedTime == 0) {
                chunk[i] = 0;
                continue;
            }
// tmpp+=expectedTime;
            v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
            chunk[i] = expectedTime * k * k;
            chunk[i] += outDegree * outDegree + outDegree*k*outDegree;
            sumCompute += chunk[i];
            sampleTotalTimes += expectedTime;
        }
// printf("expected:%llu\n", tmpp);
        e_size chunkSz = sumCompute / threads;

        v_size p = 0;
        e_size tmp = 0;
        for(v_size i = 0; i < sz; i++) {
            tmp += chunk[i];
            if(tmp >= p * chunkSz) {
                pChunk[p] = i;
                p++;
                if(p == threads) break;
            }
        }
        while(p <= threads) {
            pChunk[p] = sz;
            p++;
        }

        v_size * pt = new v_size[threads];
        for(v_size i = 0; i < threads; i++) {
            pt[i] = pChunk[i];
            printf("%u %u\n", i, pt[i]);
        }

        double ans = 0.0;

        #pragma omp parallel reduction(+:t, ans)
        {
            int threadId = omp_get_thread_num();
            std::uniform_real_distribution<double> uiDistribution(0, 1);
            v_size i = pChunk[threadId];

            while((i = __sync_fetch_and_add(&pt[threadId], 1)) < pChunk[threadId + 1]) {
                v_size u = nodes[i];

                e_size expectedSampleTime
                    = std::round(sampleTimes * (experiments[i] / sumW) + 0.000001);

                if(expectedSampleTime == 0) continue;

                sortByColor[threadId] = g->pEdge + g->pIdx2[u];
                // sortGraph(u);
                double sumD = computeDP(u, threadId);

                v_size tt = 0;
                for(v_size j = 0; j < expectedSampleTime; j++) {
                    tt += sampleOneTime(i, u, uiDistribution, e[threadId], threadId);
                }
                t += tt;
                ans += 1.0*tt/expectedSampleTime*experiments[i];
                // sampleTotalTimes += expectedSampleTime;
            }

            for(int j = 0; j < threads; j++) {
                if(j != threadId && pt[j] < pChunk[j + 1]) {
                    v_size i;
                    while((i = __sync_fetch_and_add(&pt[j], 1)) < pChunk[j + 1]) {
                        v_size u = nodes[i];
                        e_size expectedSampleTime 
                            = std::round(sampleTimes * experiments[i] / sumW);
                        if(expectedSampleTime == 0) {
                            continue; 
                        }
                        sortByColor[threadId] = g->pEdge + g->pIdx2[u];
                        computeDP(u, threadId);
                        v_size tt = 0;
                        for(v_size l = 0; l < expectedSampleTime; l++) {
                            tt += sampleOneTime(i, u, uiDistribution, e[threadId], threadId);
                        }
                        t += tt;
                        // __sync_fetch_and_add(&sampleTotalTimes, expectedSampleTime);
                    }
                }
            }
        }
        
        if(sampleTotalTimes == 0) {
            std::discrete_distribution<int> 
              udistribution(experiments, experiments + sz);
            std::uniform_real_distribution<double> uiDistribution(0, 1);
printf("|not expected %llu \n", sampleTimes - sampleTotalTimes);
            while(sampleTotalTimes < sampleTimes) {
                int id = udistribution(e[0]);
                v_size u = nodes[id];
                sortByColor[0] = g->pEdge + g->pIdx2[u];
                // sortGraph(u);
                computeDP(u, 0);
                t += sampleOneTime(id, u, uiDistribution, e[0], 0);
                sampleTotalTimes++;
            }
             // printf("|small %.6f %u %u", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
             // return 1.0 * t / sampleTotalTimes * sumW;
        }
        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
// printf("nz : %d\n", nz);
        printf("| %.6f %u %llu\n", 1.0 * t / sampleTotalTimes, t, sampleTotalTimes);
        // printf("| %.8f", expectedN / sumW);
        delete [] chunk;
        delete [] pChunk;
        delete [] pt;

        return 1.0 * t / sampleTotalTimes * sumW;
        //return ans;
    }
};

#endif