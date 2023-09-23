#ifndef PIVOTERX_HPP
#define PIVOTERX_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"
#include <cassert>
#include <assert.h>
#include <tuple>
#include <vector>
#include <chrono> 
using std::tuple;
using Pair = std::pair<v_size, v_size>;
using std::vector;

#define pP first
#define pR second

class PivoterX {
private:
    Graph * g;
    v_size k;
    c_size ** cnt = nullptr;
    c_size ** C = nullptr;
    v_size ** tmp;
    hopstotchHash * hashTable;
    v_size nodeNum = 0;

    std::vector<std::vector<bool>> isClique;

    v_size branches = 0;

public:
    PivoterX(Graph * g_) { g = g_;}
    ~PivoterX() {
        for(v_size i = 0; i <= g->degeneracy; i++) {
            delete [] C[i];
            delete [] cnt[i];
            delete [] tmp[i];
        }
        delete [] tmp;
        delete [] cnt;
        delete [] C;
        delete [] hashTable;
    }
// #define DEBUG
    void runV(v_size k_, v_size deb=111111111) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
        hashTable = new hopstotchHash[g->vCnt];
        for(v_size u = 0; u < g->vCnt; u++) {
            // printf("%u\n", u);fflush(stdout);
            if(g->pIdx[u + 1] == g->pIdx[u]) continue;
            hashTable[u].build(g->pEdge + g->pIdx[u], g->pIdx[u + 1] - g->pIdx[u]);
        }
    // printf("build hash tables\n");fflush(stdout);
#ifdef DEBUG
g->print();
#endif

        k = k_;
        cnt = new c_size*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++)
            cnt[i] = new c_size[k + 1]();
        tmp = new v_size*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++)
            tmp[i] = new v_size[g->degeneracy]();
        C = new c_size*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++) {
            C[i] = new c_size[k + 1]();
        }
        C[0][0] = 1.0;
        C[1][0] = 1.0;
        C[1][1] = 1.0;
        for(v_size i = 2; i <= g->degeneracy; i++) {
            C[i][0] = 1.0;
            if(i < k + 1) C[i][i] = 1.0;
            for(v_size j = 1; j < i && j < k + 1; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
        LinearSet * S = new LinearSet(g, hashTable);
        isClique.resize(g->degeneracy + 1);
        for(v_size i = 1; i <= g->degeneracy; i++) {
            isClique[i].resize(g->vCnt);
        }


        double timeStape = clock();
        std::fill(cnt[0], cnt[0] + k+1, 0.0);
        cnt[0][1] = g->vCnt;
        
        for(v_size u = 0; u < g->vCnt; u++) {
// double timeStape = clock();
// auto start = std::chrono::high_resolution_clock::now();
            Pair section = S->sort(u, {0, g->vCnt});
            search3(section, S, 1);
            // searchWithoutPivot(section, S, 1);
            for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
                cnt[0][i] += cnt[1][i-1];
            }
// auto finish = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> elapsed = finish - start;
// printf("%u:%.10f s\n", u, elapsed.count() );
// printf("%u:%.10f s\n", u, (clock() - timeStape) / CLOCKS_PER_SEC );
        }

        printf("%f s\n", (clock() - timeStape) / CLOCKS_PER_SEC );
        for(v_size i = 2; i <= k; i++)
            printf("%d-clique: %.2f\n", i, cnt[0][i]);

        printf("branches:%u\n", branches);
        fflush(stdout);

    }

    void search3(Pair section, LinearSet * S, v_size d = 0) {
branches++;
#ifdef DEBUG
printf("deep %u:", d);
for(v_size i = 0; i < section.pR; i++) {
    printf("%u ", (*S)[i]);
}printf("\n");
#endif
        std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
        cnt[d][1] = section.pR - section.pP;
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
            cnt[d][2] = c_size(f1 + f2 + f3);
            if(f1 && f2 && f3) cnt[d][3] = 1.0;
            return;
        }

        v_size pivot, pivotDeg=0, pivotDeg2, cliqueSize = 0, partitalCliqueSize = 0;
        //find the pivot in C and X, update tmp[d], pivot and pivotDeg
        bool isCC = S->findPivotsAndCliqueX3(section, tmp[d], 
            pivot, pivotDeg, pivotDeg2, cliqueSize, partitalCliqueSize);
        //A=pivotDeg-partitalCliqueSize=pivotDeg2; 
        //B=sec.pR-cliqueSize-pivotDeg2; 
        //C=partitalCliqueSize; 
        //D=cliqueSize-partitalCliqueSize
#ifdef DEBUG
printf("pivot %u, pd %u, pd2 %u\n", pivot, pivotDeg, pivotDeg2);
printf("cs %u, pcs %u\n", cliqueSize, partitalCliqueSize);
#endif
        if(isCC) {
            for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
                cnt[d][i] = C[pivotDeg + 1][i];
            }
            return;
        }
        
        //A+C
        search3({section.pP, section.pP + pivotDeg}, S, d+1);
        for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
            cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
        }
#ifdef DEBUG
printf("back to deep %u:\n", d);
printf("pivot cnt[d][3] %.0f\n", cnt[d][3]);
#endif
        if(cliqueSize > 0) {//C+D
            for(v_size i = 2; i <= k && i <= cliqueSize; i++) {
                cnt[d][i] += C[cliqueSize][i];
            }
#ifdef DEBUG
printf("clique add cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            if(partitalCliqueSize > 0) {//-C
                for(v_size i = 2; i <= partitalCliqueSize && i <= k; i++) {
                    cnt[d][i] -= C[partitalCliqueSize][i];
                }
#ifdef DEBUG
printf("partitalCliqueSize minus cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            }
        }

        //enumerate B
        v_size ed = section.pR - cliqueSize - pivotDeg2;
        for(v_size i = 0; i < ed; i++) {
            v_size v = tmp[d][i];
#ifdef DEBUG
printf("Bdeep %u, i%u, h %u\n", d, i, v);
#endif
            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            search3(sec, S, d+1);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
#ifdef DEBUG
printf("back to Bdeep %u:\n", d);
printf("cnt[d][3] %.0f\n", cnt[d][3]);
#endif
        }

        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++)
            isClique[d][tmp[d][i]] = true;
        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++) {
            v_size v = tmp[d][i];
            if(hashTable[pivot].contain(v)) continue;
            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            v_size cSize = 0, aSize = 0;
            for(v_size j = sec.pP; j < sec.pR; j++) {
                if(!hashTable[pivot].contain((*S)[j])) cSize++;
                else if(!isClique[d][(*S)[j]]) 
                    tmp[d][ed2 + aSize++] = (*S)[j];
            }
            if(aSize == 0) continue;
#ifdef DEBUG
printf("Ddeep %u, i%u, h %u A:", d, i, v);
for(v_size j = 0; j < aSize; j++) {
    v_size w = tmp[d][j + ed2];
    printf("%u ", w);
}printf("\n");
#endif

            for(v_size j = 0; j < aSize; j++) {
                v_size w = tmp[d][j + ed2];
#ifdef DEBUG
printf("Ddeep %u, i%u, j %u h %u, w %u\n", d, i, j, v, w);
#endif
                S->changeTo(w, --sec.pR);
                v_size newPR = 0;
                for(v_size j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[w].contain((*S)[j])) {
                        S->changeToByPos(j, newPR++);
                    }
                }
                search3({sec.pP, newPR}, S, d+1);
                cnt[d][2] += 1;
                for(v_size i = 3; i <= k && i-2 <= newPR - sec.pP; i++) {
                    cnt[d][i] += cnt[d+1][i-2];
                    if(cnt[d+1][i-2] == 0) break;
                }
#ifdef DEBUG
printf("back to Ddeep %u:\n", d);
printf("cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            }
        }
        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++)
            isClique[d][tmp[d][i]] = false;
    }

    void search2(Pair section, LinearSet * S, v_size d = 0) {
branches++;
#ifdef DEBUG
printf("deep %u:", d);
for(v_size i = 0; i < section.pR; i++) {
    printf("%u ", (*S)[i]);
}printf("\n");
#endif
        std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
        cnt[d][1] = section.pR - section.pP;
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
            cnt[d][2] = c_size(f1 + f2 + f3);
            if(f1 && f2 && f3) cnt[d][3] = 1.0;
            return;
        }

        v_size pivot, pivotDeg=0, pivotDeg2, cliqueSize = 0, partitalCliqueSize = 0;
        //find the pivot in C and X, update tmp[d], pivot and pivotDeg
        bool isCC = S->findPivotsAndCliqueX(section, tmp[d], 
            pivot, pivotDeg, pivotDeg2, cliqueSize, partitalCliqueSize);
        //A=pivotDeg-partitalCliqueSize=pivotDeg2; 
        //B=sec.pR-cliqueSize-pivotDeg2; 
        //C=partitalCliqueSize; 
        //D=cliqueSize-partitalCliqueSize
#ifdef DEBUG
printf("pivot %u, pd %u, pd2 %u\n", pivot, pivotDeg, pivotDeg2);
printf("cs %u, pcs %u\n", cliqueSize, partitalCliqueSize);
#endif
        if(isCC) {
            for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
                cnt[d][i] = C[pivotDeg + 1][i];
            }
            return;
        }
        
        //A+C
        search2({section.pP, section.pP + pivotDeg}, S, d+1);
        for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
            cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
        }
#ifdef DEBUG
printf("back to deep %u:\n", d);
printf("pivot cnt[d][3] %.0f\n", cnt[d][3]);
#endif
        if(cliqueSize > 0) {//C+D
            for(v_size i = 2; i <= k && i <= cliqueSize; i++) {
                cnt[d][i] += C[cliqueSize][i];
            }
#ifdef DEBUG
printf("clique add cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            if(partitalCliqueSize > 0) {//-C
                for(v_size i = 2; i <= partitalCliqueSize && i <= k; i++) {
                    cnt[d][i] -= C[partitalCliqueSize][i];
                }
#ifdef DEBUG
printf("partitalCliqueSize minus cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            }
        }

        //enumerate B
        v_size ed = section.pR - cliqueSize - pivotDeg2;
        for(v_size i = 0; i < ed; i++) {
            v_size v = tmp[d][i];
#ifdef DEBUG
printf("Bdeep %u, i%u, h %u\n", d, i, v);
#endif
            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            search2(sec, S, d+1);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
#ifdef DEBUG
printf("back to Bdeep %u:\n", d);
printf("cnt[d][3] %.0f\n", cnt[d][3]);
#endif
        }
//         for(v_size i = ed; i < ed + cliqueSize - partitalCliqueSize; i++) {
//             v_size v = tmp[d][i];
//             section.pR--;
//             S->changeTo(v, section.pR);
//             Pair sec = S->sort2(v, section);

//             v_size cSize = 0;
//             for(v_size j = sec.pP; j < sec.pR; j++) {
//                 if(!hashTable[pivot].contain((*S)[j])) {
//                     S->changeToByPos(j, cSize);
//                     cSize++;
//                 }
//             }
// for(v_size j = 0; j < cSize; j++) {
//     for(v_size l = j + 1; l < cSize; l++) {
//         assert(hashTable[(*S)[j]].contain((*S)[l]));
//     }
// }
// for(v_size j = cSize; j < sec.pR; j++) {
//     bool f = true;
//     for(v_size l = 0; l < cSize; l++) {
//         if(!hashTable[(*S)[l]].contain((*S)[j])) {
//             f = false;
//             break;
//         }
//     }
    
//     assert(!f || hashTable[pivot].contain((*S)[j]));
// }
//             if(cSize+sec.pP == sec.pR) continue;

//             search2(sec, S, d+1);
//             for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
//                 cnt[d][i] += cnt[d+1][i-1];
//             }
//             for(v_size i = 2; i <= k && i <= cSize+1; i++) {
//                 cnt[d][i] -= C[cSize][i-1];
//             }
//         }
        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++)
            isClique[d][tmp[d][i]] = true;
        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++) {
            v_size v = tmp[d][i];
            if(hashTable[pivot].contain(v)) continue;
            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            v_size cSize = 0, aSize = 0;
            for(v_size j = sec.pP; j < sec.pR; j++) {
                if(!hashTable[pivot].contain((*S)[j])) cSize++;
                else if(!isClique[d][(*S)[j]]) 
                    tmp[d][ed2 + aSize++] = (*S)[j];
            }
            if(aSize == 0) continue;
#ifdef DEBUG
printf("Ddeep %u, i%u, h %u A:", d, i, v);
for(v_size j = 0; j < aSize; j++) {
    v_size w = tmp[d][j + ed2];
    printf("%u ", w);
}printf("\n");
#endif

            for(v_size j = 0; j < aSize; j++) {
                v_size w = tmp[d][j + ed2];
#ifdef DEBUG
printf("Ddeep %u, i%u, j %u h %u, w %u\n", d, i, j, v, w);
#endif
                S->changeTo(w, --sec.pR);
                v_size newPR = 0;
                for(v_size j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[w].contain((*S)[j])) {
                        S->changeToByPos(j, newPR++);
                    }
                }
                search2({sec.pP, newPR}, S, d+1);
                cnt[d][2] += 1;
                for(v_size i = 3; i <= k && i-2 <= newPR - sec.pP; i++) {
                    cnt[d][i] += cnt[d+1][i-2];
                    if(cnt[d+1][i-2] == 0) break;
                }
#ifdef DEBUG
printf("back to Ddeep %u:\n", d);
printf("cnt[d][3] %.0f\n", cnt[d][3]);
#endif
            }
        }
        for(v_size i = ed, ed2 = ed + cliqueSize; i < ed2; i++)
            isClique[d][tmp[d][i]] = false;
    }
    
    void search(Pair section, LinearSet * S, v_size d = 0) {
        std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
        cnt[d][1] = section.pR - section.pP;
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
            cnt[d][2] = c_size(f1 + f2 + f3);
            if(f1 && f2 && f3) cnt[d][3] = 1.0;
            return;
        }

        v_size pivot, pivotDeg=0, pivotDeg2, cliqueSize = 0, partitalCliqueSize = 0;
        //find the pivot in C and X, update tmp[d], pivot and pivotDeg
        bool isCC = S->findPivotsAndCliqueX(section, tmp[d], 
            pivot, pivotDeg, pivotDeg2, cliqueSize, partitalCliqueSize);
        //A=pivotDeg-partitalCliqueSize=pivotDeg2; 
        //B=sec.pR-cliqueSize-pivotDeg2; 
        //C=partitalCliqueSize; 
        //D=cliqueSize-partitalCliqueSize

        if(isCC) {
            for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
                cnt[d][i] = C[pivotDeg + 1][i];
            }
            return;
        }
        
        //A+C
        search({section.pP, section.pP + pivotDeg}, S, d+1);
        for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
            cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
        }

        //enumerate B
        v_size ed = section.pR - cliqueSize - pivotDeg2;
        for(v_size i = 0; i < ed; i++) {
            v_size v = tmp[d][i];

            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            search(sec, S, d+1);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
        }
        for(v_size i = ed; i < ed + cliqueSize - partitalCliqueSize; i++) {
            v_size v = tmp[d][i];
            section.pR--;
            S->changeTo(v, section.pR);
            Pair sec = S->sort2(v, section);

            search(sec, S, d+1);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
        }
    }
    

};

#undef pP
#undef pR

#endif
