
#ifndef LINEARSET
#define LINEARSET
#include <tuple>
#include <utility>
#include "type.hpp"
#include "hopstotchHash.hpp"
#include "../graph/graph.hpp"

using Pair = std::pair<v_size, v_size>;

#define pP first
#define pR second


class LinearSet {
private:
    v_size * vSet;
    v_size * fIndex;
    Graph * g;
    hopstotchHash * hashTable;
    bool * vis;
    
public:
    LinearSet(Graph * g_, hopstotchHash * hashTable_) {
        g = g_;
        hashTable = hashTable_;
        vSet = new v_size[g->vCnt];
        fIndex = new v_size[g->vCnt];
        for(v_size i = 0; i < g->vCnt; i++) {
            vSet[i] = fIndex[i] = i;
        }
        vis = new bool[g->degeneracy];
    }

    ~LinearSet() { delete [] fIndex; delete [] vSet; delete [] vis;}
    v_size * begin() {
        return vSet;
    }
    bool isIn(v_size v, v_size l, v_size r) {
        return l <= fIndex[v] && fIndex[v] < r;
    }

    v_size operator [] (v_size i) {
        // if(i >= g->vCnt) {
        //     printf("error index\n"); return -1;
        // }
        return vSet[i];
    }

    void changeTo(v_size u, v_size p) {
        v_size pU = fIndex[u];
        std::swap(fIndex[u], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    void changeToByPos(v_size pU, v_size p) {
        std::swap(fIndex[vSet[pU]], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    Pair sort(v_size u, const Pair & sec) {
        v_size l, r;
        l = r = sec.pP;

        for(v_size i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
            v_size v = g->pEdge[i];
            // if(sec.pP <= fIndex[v] && fIndex[v] < sec.pR) {
            changeTo(v, r++);
            // }
        }

        return {l, r};
    }

    Pair sort2(v_size u, const Pair & sec) {
        v_size i = sec.pP, j = sec.pP;
        while(i < sec.pR && j < sec.pR) {
            if(hashTable[u].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        return {sec.pP, j};
    }

    Pair sort3(v_size u, const Pair & sec, v_size p, v_size & newP_) {
        //clique
        v_size l = p;
        for(v_size i = p - 1; i >= sec.pP; i--) {
            if(hashTable[u].contain(vSet[i])) {
                changeToByPos(i, --l);
            }
        }

        v_size i = p, j = p;
        while(i < sec.pR) {
            if(hashTable[u].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        v_size newP = p;
        for(v_size i = p; i < j; i++) {
            bool f = true;
            for(v_size j = l; j < newP; j++) {
                if(!hashTable[vSet[i]].contain(vSet[j])) {
                    f = false; break;
                }
            }
            if(f) {
                changeToByPos(i, newP++);
            }
        }
        newP_ = newP;
        return {l, j};
    }

    v_size findClique(const Pair & sec) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;
                }
            }
        }

        changeTo(pivot, sec.pP);

        v_size l = sec.pP, r = sec.pR;
        do {
            for(v_size i = l + 1; i < r; ) {
                if(!hashTable[vSet[l]].contain(vSet[i])) {
                    changeToByPos(i, --r);
                }
                else i++;
            }
            
            l++;
        }while(l < r);
        return r;
    }

    Pair mtMaxClique(const Pair & sec, v_size p) {
        //p >= 2
        v_size pivot = vSet[sec.pP];
        // changeToByPos(p - 1, sec.pR - 1);

        v_size j = p;
        for(v_size i = p; i < sec.pR; i++) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeToByPos(i, j++);
            }
        }
        
        v_size newP = p;
        if(p == sec.pP + 1) newP += 1;
        for(v_size i = newP; i < j; i++) {
            bool f = true;
            for(v_size j = sec.pP + 1; j < newP; j++) {
                if(!hashTable[vSet[i]].contain(vSet[j])) {
                    f = false; break;
                }
            }
            if(f) changeToByPos(i, newP++);
        }

        return {newP, j};
    }

    void findPivot(const Pair & sec, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;
                }
            }
        }

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }
// printf("%u %u %u %u\n", j, i, sec.pR, vSet[sec.pR - 1]);
        pivot_ = pivot;
        pivotDeg_ = pivotDeg;
// if(sec.pP > g->vCnt) {
//     printf("%u\n", sec.pP);fflush(stdout);
// }
        // v_size pivot = vSet[sec.pP];
        // v_size pivotDeg = sec.pR - sec.pP;
        // v_size mid = g->degeneracy;
        // memset(vis, false, sizeof(bool)*(sec.pR - sec.pP));

        // for(auto i = sec.pP; i < sec.pR; i++) {
        //     v_size v = vSet[i];
        //     v_size tmp = sec.pP, tmp2 = 0;

        //     for(auto j = sec.pP; j < sec.pR; j++) {
        //         v_size u = vSet[j];
        //         if(hashTable[u].contain(v)) {
        //             // changeToByPos(j, tmp++);
        //             tmp++;
        //             for(auto k = j + 1; k < sec.pR; k++)  {
        //                 v_size w = vSet[k];
        //                 if(!vis[k - sec.pP] && hashTable[u].contain(w)) {
        //                     vis[k - sec.pP] = true;
        //                     tmp2++;
        //                 }
        //             }
        //         }
        //     }

        //     if(tmp2 < mid) {
        //         mid = tmp2;
        //         pivot = v;
        //         pivotDeg = tmp;
        //     }
        // }

        // changeTo(pivot, sec.pR - 1);

        // v_size i = sec.pP, j = sec.pP;
        // while(j < sec.pP + pivotDeg) {
        //     if(hashTable[pivot].contain(vSet[i])) {
        //         changeTo(vSet[i], j++);
        //     }
        //     i++;
        // }

        // pivot_ = pivot;
        // pivotDeg_ = pivotDeg;
    }
    
    int findPivotAndCopyToTmpMem(const Pair & sec, v_size * tmpP, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;
        int num = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp; num = 0;
                }
                else if(tmp == pivotDeg){
                    num++;
                }
            }
        }

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg;

        memcpy(tmpP, vSet + sec.pP, sizeof(v_size) * (sec.pR - sec.pP));

        return num;
    }

    v_size findDensity(const Pair & section) {
        v_size n = 0;
        for(v_size i = section.pP; i < section.pR; i++) {
            for(v_size j = i + 1; j < section.pR; j++) {
                if(hashTable[vSet[i]].contain(vSet[j])) {
                    n++;
                }
            } 
        }
        return 2*n/section.pR;
    }

    int findPivotsAndClique(const Pair & sec, v_size * tmpP, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;
        int num = 0;

        // int a[20];
        // bool b[20];
        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp; num = 0;
                }
                else if(tmp == pivotDeg){
                    // a[num++] = v;
                    num++;
                }
            }
        }

        // for(v_size i = 0; i < num; i++) {
        //     if(hashTable[pivot].contain(a[i])) {
        //         b[i] = true;
        //     }
        //     else b[i] = false;
        // }
        // for(v_size i = 0; i < num; i++) {
        //     if(!b[i]) continue;
        //     for(v_size j = i + 1; j < num; j++) {
        //         if(!b[i]) continue;
        //         if(!hashTable[a[i]].contain(a[j])) b[j] = false;
        //     }
        // }
        // int cnt = 0;
        // for(v_size i = 0; i < num; i++) {
        //     if(b[i]) cnt++;
        // }
        // printf("c %d ", cnt);

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg;

        memcpy(tmpP, vSet + sec.pP, sizeof(v_size) * (sec.pR - sec.pP));

        return num;
    }

//     bool findPivotsAndCliqueX(Pair & sec, v_size * tmpP, v_size & pivot_, v_size & pivotDeg_, v_size & cliqueSize) {
//         //find a large clique at first

        
//         v_size pivot = vSet[sec.pP], pp = sec.pP;
//         v_size pivotDeg = 0;
//         int num = 0;

//         // int a[20];
//         // bool b[20];
//         for(auto i = sec.pP; i < sec.pR; i++) tmpP[i] = 0;
//         for(auto i = sec.pP; i < sec.pR; i++) {
//             v_size v = vSet[i];
//             if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
//                 v_size tmp = 0;

//                 for(auto j = sec.pP; j < sec.pR; j++) {
//                     if(hashTable[v].contain(vSet[j]) ) tmp++;
//                 }
                
//                 tmpP[i] = tmp;
//                 if(tmp > pivotDeg) {
//                     pivot = v; pivotDeg = tmp; num = 0;
//                     pp = i;
//                 }
//                 else if(tmp == pivotDeg){
//                     // a[num++] = v;
//                     num++;
//                 }
//             }
//         }

//         changeTo(pivot, --sec.pR);
//         std::swap(tmpP[pp], tmpP[sec.pR]);

//         v_size i = sec.pP, j = sec.pP;
//         while(j < sec.pP + pivotDeg) {
//             if(hashTable[pivot].contain(vSet[i])) {
//                 std::swap(tmpP[i], tmpP[j]);
//                 changeTo(vSet[i], j++);
//             }
//             i++;
//         }

//         pivot_ = pivot;
//         pivotDeg_ = pivotDeg;

//         //find a large clique
//         // j = sec.pP + pivotDeg;
//         if(sec.pR > sec.pP + pivotDeg) {
//             for(auto i = sec.pP + pivotDeg + 1; i < sec.pR; i++) {
//                 if(tmpP[i] > tmpP[j]) {
//                     j = i;
//                 }
//             }
//             cliqueSize++;
//             v_size u = vSet[j];
//             changeToByPos(j, sec.pR - 1);
//             std::swap(tmpP[sec.pR - 1], tmpP[j]);

//             j = sec.pP + pivotDeg;
//             for(auto i = sec.pP + pivotDeg; i < sec.pR - 1; i++) {
//                 if(hashTable[u].contain(vSet[i]) ) {
//                     std::swap(tmpP[i], tmpP[j]);
//                     changeToByPos(i, j++);
//                 }
//             }

//             while(j > sec.pP + pivotDeg) {
//                 v_size pMax = sec.pP + pivotDeg;
//                 for(auto i = sec.pP + pivotDeg + 1; i < j; i++) {
//                     if(tmpP[i] > tmpP[pMax]) pMax = i;
//                 }

//                 cliqueSize++;
//                 u = vSet[pMax];
//                 changeToByPos(pMax, --j);
//                 std::swap(tmpP[j], tmpP[pMax]);

//                 changeToByPos(j, sec.pR - cliqueSize);
//                 // std::swap(tmpP[j], tmpP[sec.pR]);

//                 v_size newJ = sec.pP + pivotDeg;
//                 for(auto i = sec.pP + pivotDeg; i < j; i++) {
//                     if(hashTable[u].contain(vSet[i])) {
//                         std::swap(tmpP[i], tmpP[newJ]);
//                         changeToByPos(i, newJ++);
//                     }
//                 }
//                 j = newJ;
//             }
// //check if find a clique
// // for(i = sec.pR - cliqueSize; i < sec.pR; i++) {
// //     for(j = i + 1; j < sec.pR; j++) {
// //         if(!hashTable[vSet[i]].contain(vSet[j])) {
// //             printf("error, not a large clique\n");fflush(stdout);
// //         }
// //     }
// // }
//         }

//         memcpy(tmpP, vSet + sec.pP, sizeof(v_size) * (sec.pR - sec.pP));

//         return num == pivotDeg && sec.pR == pivotDeg;
//     }
// #define DEBUG
    bool findPivotsAndCliqueX(Pair & sec, v_size * tmpP, 
        v_size & pivot_, v_size & pivotDeg_, v_size & pivotDeg2,
        v_size & cliqueSize, v_size & partitalCliqueSize) {
        v_size pivot = vSet[sec.pP], pp = sec.pP;
        v_size pivotDeg = 0, num = 0;

        for(auto i = sec.pP; i < sec.pR; i++) tmpP[i] = 0;
        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }
                
                tmpP[i] = tmp;
                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;// num = 0;
                    pp = i;
                }
                else if(tmp == pivotDeg){
                    num++;
                }
            }
        }
        
        if((v_size)num == pivotDeg && sec.pR == pivotDeg + 1) {
            pivot_ = pivot;
            pivotDeg_ = pivotDeg;
            return true;
        }

        cliqueSize++;
        changeTo(pivot, sec.pR - cliqueSize);
        std::swap(tmpP[pp], tmpP[sec.pR - cliqueSize]);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                std::swap(tmpP[i], tmpP[j]);
                changeTo(vSet[i], j++);
            }
            i++;
        }

        while(j > 0) {
            v_size pMax = 0;
            for(v_size i = 1; i < j; i++) {
                if(tmpP[i] > tmpP[pMax]) pMax = i;
            }

            cliqueSize++;
            pivot = vSet[pMax];
            changeToByPos(pMax, --j);
            std::swap(tmpP[j], tmpP[pMax]);

            changeToByPos(j, sec.pR - cliqueSize);
            std::swap(tmpP[j], tmpP[sec.pR - cliqueSize]);

            v_size newJ = 0;
            for(v_size i = 0; i < j; i++) {
                if(hashTable[pivot].contain(vSet[i])) {
                    std::swap(tmpP[i], tmpP[newJ]);
                    changeToByPos(i, newJ++);
                }
            }
            j = newJ;
        }
#ifdef DEBUG
printf("find a clique:");
for(i = sec.pR - cliqueSize; i < sec.pR; i++) printf("%u ", vSet[i]);
printf("\n");
#endif
//check if find a clique
// for(i = sec.pR - cliqueSize; i < sec.pR; i++) {
//     for(j = i + 1; j < sec.pR; j++) {
//         if(!hashTable[vSet[i]].contain(vSet[j])) {
//             printf("error, not a large clique\n");fflush(stdout);
//         }
//     }
// }
// if(sec.pR <= cliqueSize) {
//     printf("pR %u, Cs %u\n", sec.pR, cliqueSize);fflush(stdout);
// }
// assert(sec.pR > cliqueSize);

        //choose a pivot apart from the clique
        pivot = vSet[sec.pP], pp = sec.pP;
        pivotDeg = 0;
        for(auto i = sec.pP; i < sec.pR - cliqueSize; i++) {
            v_size v = vSet[i];
            v_size tmp = 0;

            for(auto j = sec.pP; j < sec.pR - cliqueSize; j++) {
                if(hashTable[v].contain(vSet[j]) ) tmp++;
            }
            
            tmpP[i] = tmp;
            if(tmp > pivotDeg) {
                pivot = v; pivotDeg = tmp;
                pp = i;
            }
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg2 = pivotDeg;
// if(sec.pR < cliqueSize + 1 || sec.pR < 1) {
//     printf("sec.pR %u cliqueSize %u\n", sec.pR, cliqueSize);fflush(stdout);
// }
// assert(pp < sec.pR - cliqueSize);
// for(v_size ii = sec.pR - cliqueSize; ii < sec.pR; ii++) {
//     for(v_size jj = ii + 1; jj < sec.pR; jj++) {
//         assert(hashTable[vSet[ii]].contain(vSet[jj]));
//     }
// }
        changeToByPos(sec.pR - cliqueSize - 1, sec.pR - 1);
// for(v_size ii = sec.pR - cliqueSize-1; ii < sec.pR-1; ii++) {
//     for(v_size jj = ii + 1; jj < sec.pR-1; jj++) {
//         assert(hashTable[vSet[ii]].contain(vSet[jj]));
//     }
// }
        // std::swap(tmpP[sec.pR - cliqueSize - 1], tmpP[sec.pR - 1]);
        if(pp < sec.pR - cliqueSize - 1) {
// assert(vSet[pp] == pivot);
            changeToByPos(pp, sec.pR - 1);
// for(v_size ii = sec.pR - cliqueSize-1; ii < sec.pR-1; ii++) {
//     for(v_size jj = ii + 1; jj < sec.pR-1; jj++) {
//         assert(hashTable[vSet[ii]].contain(vSet[jj]));
//     }
// }
        }
        // std::swap(tmpP[pp], tmpP[sec.pR - 1]);
        --sec.pR;
// for(v_size ii = sec.pR - cliqueSize; ii < sec.pR; ii++) {
//     for(v_size jj = ii + 1; jj < sec.pR; jj++) {
//         assert(hashTable[vSet[ii]].contain(vSet[jj]));
//     }
// }

        i = sec.pP, j = sec.pP;
        // while(j < sec.pP + pivotDeg) {
        while(i < sec.pR - cliqueSize) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }
// assert(pivotDeg2 == j);
// for(v_size ii = 0; ii < pivotDeg2; ii++) 
//     assert(hashTable[pivot].contain(vSet[ii]));
        // memcpy(tmpP, vSet + j, 
        //     sizeof(v_size) * (sec.pR - j - cliqueSize));
        memcpy(tmpP, vSet + j, 
            sizeof(v_size) * (sec.pR - j));
// for(v_size ii = sec.pR - cliqueSize- j; ii < sec.pR- j; ii++) {
//     for(v_size jj = ii + 1; jj < sec.pR- j; jj++) {
//         assert(hashTable[tmpP[ii]].contain(tmpP[jj]));
//     }
// }
        // num = sec.pR - pivotDeg2 - cliqueSize;
        while(i < sec.pR) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
                partitalCliqueSize++;
                pivotDeg_++;
            }
            // else tmpP[num++] = vSet[i];
            i++;
        }
// assert(num == sec.pR - pivotDeg_);
// for(v_size ii = 0; ii < num; ii++) {
//     assert(!hashTable[pivot].contain(tmpP[ii]));
// }
// assert(pivotDeg_ == j);
// for(v_size ii = j; ii < sec.pR; ii++) {
//     assert(!hashTable[pivot].contain(vSet[ii]));
// }

        return false;
    }

    bool findPivotsAndCliqueX3(Pair & sec, v_size * tmpP, 
        v_size & pivot_, v_size & pivotDeg_, v_size & pivotDeg2,
        v_size & cliqueSize, v_size & partitalCliqueSize) {
        v_size pivot = vSet[sec.pP], pp = sec.pP;
        v_size pivotDeg = 0, num = 0;

        for(auto i = sec.pP; i < sec.pR; i++) tmpP[i] = 0;
        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }
                
                tmpP[i] = tmp;
                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;// num = 0;
                    pp = i;
                }
                else if(tmp == pivotDeg){
                    num++;
                }
            }
        }
        
        if((v_size)num == pivotDeg && sec.pR == pivotDeg + 1) {
            pivot_ = pivot;
            pivotDeg_ = pivotDeg;
            return true;
        }

        cliqueSize++;
        changeTo(pivot, sec.pR - cliqueSize);
        std::swap(tmpP[pp], tmpP[sec.pR - cliqueSize]);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                std::swap(tmpP[i], tmpP[j]);
                changeTo(vSet[i], j++);
            }
            i++;
        }

        while(j > 0) {
            v_size pMax = 0;
            for(v_size i = 1; i < j; i++) {
                if(tmpP[i] > tmpP[pMax]) pMax = i;
            }

            cliqueSize++;
            pivot = vSet[pMax];
            changeToByPos(pMax, --j);
            std::swap(tmpP[j], tmpP[pMax]);

            changeToByPos(j, sec.pR - cliqueSize);
            std::swap(tmpP[j], tmpP[sec.pR - cliqueSize]);

            v_size newJ = 0;
            for(v_size i = 0; i < j; i++) {
                if(hashTable[pivot].contain(vSet[i])) {
                    std::swap(tmpP[i], tmpP[newJ]);
                    changeToByPos(i, newJ++);
                }
            }
            j = newJ;
        }
#ifdef DEBUG
printf("find a clique:");
for(i = sec.pR - cliqueSize; i < sec.pR; i++) printf("%u ", vSet[i]);
printf("\n");
#endif

        //choose a pivot apart from the clique
        pivot = vSet[sec.pP], pp = sec.pP;
        pivotDeg = 0;
        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            v_size tmp = 0;

            for(auto j = sec.pP; j < sec.pR - cliqueSize; j++) {
                if(hashTable[v].contain(vSet[j]) ) tmp++;
            }
            
            tmpP[i] = tmp;
            if(tmp > pivotDeg) {
                pivot = v; pivotDeg = tmp;
                pp = i;
            }
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg2 = pivotDeg;

        changeToByPos(sec.pR - cliqueSize - 1, sec.pR - 1);

        // std::swap(tmpP[sec.pR - cliqueSize - 1], tmpP[sec.pR - 1]);
        if(pp < sec.pR - cliqueSize - 1) {
// assert(vSet[pp] == pivot);
            changeToByPos(pp, sec.pR - 1);
        }
        // std::swap(tmpP[pp], tmpP[sec.pR - 1]);
        --sec.pR;

        i = sec.pP, j = sec.pP;
        // while(j < sec.pP + pivotDeg) {
        while(i < sec.pR - cliqueSize) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }
        memcpy(tmpP, vSet + j, 
            sizeof(v_size) * (sec.pR - j));
        while(i < sec.pR) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
                partitalCliqueSize++;
                pivotDeg_++;
            }
            // else tmpP[num++] = vSet[i];
            i++;
        }

        return false;
    }

    void copy(v_size * tmpP, const Pair & sec) {
        for(v_size i = 0; i < sec.pR - sec.pP; i++) {
            v_size v = tmpP[i];
            changeTo(v, i);
        }
    }
};

#endif