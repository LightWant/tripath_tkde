#ifndef MULTILEVEL_HPP
#define MULTILEVEL_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/bitmap.hpp"
#include "../tools/linearSet.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;
#define pP first
#define pR second

class multiLevel {
private:
    Graph * g, *subG;
    v_size k;
    double ** cnt;
    c_size ** C = nullptr;
    hopstotchHash * hashTable;
    v_size len = 1000000;
    v_size edges;
    LinearSet * S;
    v_size ** tmp;
    v_size maxDeepth = 0;

    //multiRun:
    v_size maxK;
    v_size * eCnt = nullptr;
    v_size * vCnt = nullptr;
    double ** ansTmp = nullptr;
    double * ansTmpMemoryPool = nullptr;

public:
    multiLevel(Graph * g_) {g = g_;}
    ~multiLevel() {
        for(v_size i = 0; i <= g->degeneracy; i++) {
            delete [] C[i];
            // delete [] tmp[i];
        }
        delete S;
        delete [] tmp[0];
        delete [] tmp;
        delete [] C;
        delete [] cnt[0];

        delete [] cnt;
        delete [] hashTable;
        delete subG;
        if(eCnt != nullptr) delete [] eCnt;
        if(vCnt != nullptr) delete [] vCnt;
        if(ansTmp != nullptr) delete [] ansTmp;
        if(ansTmpMemoryPool != nullptr) delete [] ansTmpMemoryPool;
    }
    void previousWork();
    void buildSubGraphFalse(v_size u);
    void buildSubGraphTrue(v_size u, Bitmap*subGraph, v_size &l, 
        e_size & eCnt, v_size & vCnt, double & experiments);
    bool cc(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    void shadowFinder(v_size l);
    void search(Pair section, LinearSet * S, v_size d = 0, v_size h=1);
    void searchForSpecificK(v_size h, v_size p, Pair section, LinearSet * S, v_size k, v_size d = 0);

    void testMultiDensity(v_size k_, v_size deb, double exactN, double a, e_size N=50000000);
};


void multiLevel::testMultiDensity(v_size k_, v_size deb, double exCnt, double alpha, e_size N) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
    printf("tripath\n");
    printf("k %u\nalpha %.1f\nN %llu\n", k, alpha, N);
    previousWork();
    g->colorG();
    
    double t = clock(), t1, tS;
    tS = t;

    v_size st = 0, ed = g->vCnt;
    
    tripath * stpObj = new tripath();
    stpObj->initForSingleNode(k-1, g, hashTable);

    std::fill(cnt[0], cnt[0] + k+1, 0.0);
    cnt[0][1] = g->vCnt;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    for(v_size u = st; u < ed; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        buildSubGraphFalse(u);
        
        if(subG->eCnt / 2 > alpha* (k-1) * subG->vCnt) {
            if(sz == nodes.size()) {
                nodes.resize(sz + ed - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {
            Pair section = S->sort(u, {0, g->vCnt});
            cnt[1][k-1] = 0;
            search(section, S, 1, 1);
            cnt[0][k] += cnt[1][k-1];

            // Pair section = S->sort(u, {0, g->vCnt});
            // searchForSpecificK(1, 0, section, S, k, 0);
        }
    }

    t1 = clock();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    printf("timesExact %.2f\n", (t1 - t)/CLOCKS_PER_SEC);
    t = t1;

    printf("dense size %u\n", sz);
    printf("maxDeepthOfPivoter %u\n", maxDeepth);
    
    double exactCnt = cnt[0][k];

    printf("exact part count %.1f\n", exactCnt);
    printf("total count %.1f\n", exCnt);
   
    if(sz > 0) {//tripath
        stpObj->init(sz, nodes);

        double sampleCnt_ = stpObj->sample(nodes, N, exCnt - exactCnt);
        printf("count of path %.6f\n", stpObj->sumW);
        printf("clique density %.6f\n", ((exCnt - exactCnt) / stpObj->sumW));
        printf("expected sparse part %.0f\n", sampleCnt_);
        printf("expected sparse part clique density %.0f\n", (sampleCnt_ / stpObj->sumW));
        double t1 = clock();
        // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
        printf("exact/total %.8f%%\n", ((exCnt - exactCnt) / exCnt)*100);
        double totalCnt = exactCnt + sampleCnt_;
        
        printf("expected answer %.0f\n", totalCnt);
        printf("exact/answer %.2f%%\n", (exactCnt / totalCnt)*100);
        // double totalCntML = exactCnt + sampleCntMultiLayer;
        // printf("| %.2f%%", (exactCnt / totalCntML)*100);
        printf("sample time %.2f\n", (t1 - t)/CLOCKS_PER_SEC);
        printf("totalTime %.2f\n", (t1 - tS)/CLOCKS_PER_SEC);
        printf("error %.8f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
        // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// scpObj->print();
        delete stpObj;
    }

    fflush(stdout);
};


void multiLevel::previousWork() {
    hashTable = new hopstotchHash[g->vCnt];
    for(v_size u = 0; u < g->vCnt; u++) {
        if(g->pIdx[u + 1] == g->pIdx[u]) continue;
        hashTable[u].build(g->pEdge + g->pIdx[u], 
            g->pIdx[u + 1] - g->pIdx[u]);
    }

    S = new LinearSet(g, hashTable);

    cnt = new double*[g->degeneracy + 1];
    double * memoryPool2 = new double[(g->degeneracy + 1)*(k+1)];
    cnt[0] = memoryPool2;
    for(v_size i = 1; i <= g->degeneracy; i++) {
        cnt[i] = cnt[i-1] + (k + 1);
    }

    C = new c_size*[g->degeneracy + 1];
    for(v_size i = 0; i <= g->degeneracy; i++) {
        C[i] = new c_size[k + 1]();
    }
    C[0][0] = 1;
    C[1][0] = 1;
    C[1][1] = 1;
    for(v_size i = 2; i <= g->degeneracy; i++) {
        C[i][0] = 1;
        if(i < k + 1) C[i][i] = 1;
        for(v_size j = 1; j < i && j < k + 1; j++) {
            C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }

    tmp = new v_size*[g->degeneracy];
    v_size * memoryPool = new v_size[g->degeneracy*g->degeneracy];
    v_size p = 0;
    for(v_size i = 0; i < g->degeneracy; i++) {
        // tmp[i] = new v_size[g->degeneracy]();
        tmp[i] = memoryPool + p;
        p += g->degeneracy;
    }

    subG = new Graph();
    subG->pIdx = new v_size[g->degeneracy];
    subG->pEdge = new v_size[g->degeneracy*g->degeneracy];
}

void multiLevel::buildSubGraphFalse(v_size u) {
    v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

    v_size st = g->pIdx2[u], l = 0;
    for(v_size i = st; i < g->pIdx[u+1]; i++) {
        v_size v = g->pEdge[i];

        for(v_size j = st; j < g->pIdx[u+1]; j++) {
            v_size w = g->pEdge[j];
            if(cc(w, v)) l++;
        }
    }

    subG->vCnt = outDegree;
    subG->eCnt = l;
}

void multiLevel::buildSubGraphTrue(v_size u, Bitmap*subGraph, 
        v_size &l, e_size & eCnt, 
        v_size & vCnt, double & experiments) {
    v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];
    v_size edges = 0;

    v_size st = g->pIdx2[u];
    for(v_size i = st; i < g->pIdx[u+1]; i++) {
        v_size v = g->pEdge[i];
        subGraph[i-st].setNowSize(outDegree);
        subGraph[i-st].Clear();
        for(v_size j = st; j < g->pIdx[u+1]; j++) {
            v_size w = g->pEdge[j];

            if(cc(w, v)) {
                subGraph[i-st].SetBit(j - st);
                edges++;
            }
        }
    }

    vCnt = outDegree;
    eCnt = edges;
    experiments = C[outDegree][l];
}


void multiLevel::search(Pair section, LinearSet * S, v_size d, v_size h) {
    std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
    cnt[d][1] = section.pR - section.pP;
    maxDeepth = std::max(maxDeepth, d);

    if(h > k) return;

    if(section.pP == section.pR) return;
    if(section.pP + 1 == section.pR) return;
    if(section.pP + 2 == section.pR) {
        if(hashTable[(*S)[section.pP]].contain((*S)[section.pP+1])) cnt[d][2] = 1.0;
        return;
    }
    if(section.pP + 3 == section.pR) {
        int f1 = hashTable[(*S)[section.pP]].contain((*S)[section.pP+1]);
        int f2 = hashTable[(*S)[section.pP]].contain((*S)[section.pP+2]);
        int f3 = hashTable[(*S)[section.pP+1]].contain((*S)[section.pP+2]);
        cnt[d][2] = double(f1 + f2 + f3);
        if(f1 && f2 && f3) cnt[d][3] = 1.0;
        return;
    }

    v_size pivot, pivotDeg;
    int numMax = S->findPivotsAndClique(section, tmp[d], pivot, pivotDeg);

    section.pR--;

    if((v_size)numMax == pivotDeg && section.pR == pivotDeg) {
        for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
            cnt[d][i] = C[pivotDeg + 1][i];
        }
        return;
    }

    search({section.pP, section.pP + pivotDeg}, S, d+1, h);
    
    for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
        cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
    }

    v_size ed = section.pR;
    for(v_size i = section.pP + pivotDeg; i < ed; i++) {
        v_size v = tmp[d][i];
        
        section.pR--;
        S->changeTo(v, section.pR);

        Pair sec = S->sort2(v, section);
        if(sec.pR == pivotDeg && ed == section.pP + pivotDeg + 1) {
            for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
            continue;
        }

        search(sec, S, d+1, h + 1);
        for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
            cnt[d][i] += cnt[d+1][i-1];
        }
    }
};

void multiLevel::searchForSpecificK(v_size hNum, 
    v_size pNum, Pair section, LinearSet * S, v_size k, v_size d) {
    if(hNum > k) return;

    maxDeepth = std::max(maxDeepth, d);

    if(section.pP == section.pR) {
        for(v_size j = 0; j <= pNum && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum][j];
        }
        return;
    }
    if(section.pP + 1 == section.pR) {
        for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum + 1][j];
        }
        return;
    }
    if(section.pP + 2 == section.pR) {
        if(hashTable[(*S)[section.pP]].contain((*S)[section.pP+1])) {
            for(v_size j = 0; j <= pNum + 2 && hNum + j <= k; j++) {
                cnt[0][hNum + j] += C[pNum + 2][j];
            }
        }
        else {
            for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
                cnt[0][hNum + j] += C[pNum + 1][j];
            }
            for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
                cnt[0][hNum + 1 + j] += C[pNum][j];
            }
        }
        return;
    }

    v_size pivot, pivotDeg;
    // v_size * tmpP = new v_size[section.pR - section.pP];
    S->findPivotAndCopyToTmpMem(section, tmp[d], pivot, pivotDeg);
    section.pR--;
// printf("pivot, left %u,right %u, num %d\n", pivotDeg, section.pR - pivotDeg, numMax);

    if(pivotDeg == 0) {
        for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
            cnt[0][hNum + 1 + j] += C[pNum][j]*section.pR;
        }
        for(v_size j = 0; j <= pNum+1 && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum+1][j];
        }
        return;
    }

    // if(pivotDeg == (v_size)numMax && pivotDeg == section.pR) {
    //     for(v_size j = 0; j <= pNum+pivotDeg+1 && hNum + j <= k; j++) {
    //         cnt[hNum + j] += C[pNum+pivotDeg+1][j];
    //     }
    //     return;
    // }

    searchForSpecificK(hNum, pNum + 1, {section.pP, section.pP + pivotDeg}, S, k, d+1);

    // if(pivotDeg == 1) {
    //     for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
    //         cnt[hNum + 1 + j] += C[pNum][j]*(section.pR - pivotDeg - numMax);
    //     }
    //     if(numMax > 0)
    //     for(v_size j = 0; j <= pNum+1 && hNum + 1 + j <= k; j++) {
    //         cnt[hNum + 1 + j] += C[pNum+1][j]*numMax;
    //     }
    //     return;
    // }
// if(d == 0) printf("%u:", section.pR);
    v_size ed = section.pR;
    for(v_size i = pivotDeg; i < ed; i++) {
        v_size v = tmp[d][i];
        
        section.pR--;

        S->changeTo(v, section.pR);
// if(section.pR == pivotDeg) printf("1\n");
// printf("deep%u:search %u %u-%u,h %u\n", d, ed+1, hNum, pNum, section.pR);
        Pair sec = S->sort2(v, section);
        searchForSpecificK(hNum + 1, pNum, sec, S, k, d+1);
    }

}

#endif
