#ifndef KCLIST_HPP
#define KCLIST_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"
#include <cassert>
#include <tuple>
#include <vector>
using Pair = std::pair<v_size, v_size>;

class kclist
{
private:
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double ans = 0;
    LinearSet *C = nullptr;

    bool cc(v_size u, v_size v) { return hashTable[u].contain(v); }
    void listing(v_size deep, std::vector<v_size> & P, const std::vector<v_size> & C);
    void listingLS(v_size deep, v_size edC);

public:
    kclist(Graph * g, v_size k):g(g), k(k) {};
    ~kclist() {
        if(C != nullptr) delete C;
        if(hashTable != nullptr) delete [] hashTable;
    };

    void run(double exCnt);

    void runLinearset(double exCnt);
};

void kclist::runLinearset(double exCnt) {
    hashTable = new hopstotchHash[g->vCnt];
    for(v_size u = 0; u < g->vCnt; u++) {
        if(g->pIdx[u + 1] == g->pIdx[u]) continue;
        hashTable[u].build(g->pEdge + g->pIdx[u], g->pIdx[u + 1] - g->pIdx[u]);
    }

    C = new LinearSet(g, hashTable);

    for(ui u = 0; u < g->vCnt; u++) {
        ui edC = 0;
        for(ui i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
            C->changeTo(g->pEdge[i], edC++);
        }

        listingLS(1, edC);
    }

    printf("ans:%.0f\n", ans);
    printf("exCnt:%.0f\n", ans);
}
void kclist::listingLS(v_size deep, v_size edC) {
    if(deep == k - 2) {
        for(ui i = 0; i < edC; i++) {
            for(ui j = i + 1; j < edC; j++) {
                if(cc((*C)[i], (*C)[j])) ans++;
            }
        }
        ans += edC;
        return;
    }

    std::vector<v_size> cands(C->begin(), C->begin() + edC);

    for(auto v:cands) {
        ui newEdC = 0;
        C->changeTo(v, --edC);

        for(ui j = 0; j < edC; j++) {
            if(cc(v, (*C)[j])) C->changeToByPos(j, newEdC++);
        }

        listingLS(deep + 1, newEdC);
    }
}

void kclist::run(double exCnt) {
    hashTable = new hopstotchHash[g->vCnt];
    for(v_size u = 0; u < g->vCnt; u++) {
        if(g->pIdx[u + 1] == g->pIdx[u]) continue;
        hashTable[u].build(g->pEdge + g->pIdx[u], g->pIdx[u + 1] - g->pIdx[u]);
    }

    for(ui u = 0; u < g->vCnt; u++) {
        std::vector<v_size> C, P;
        for(ui i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
            C.push_back(g->pEdge[i]);
        }

        P.push_back(u);
        listing(1, P, C);
    }

    printf("ans:%.0f\n", ans);
    printf("exCnt:%.0f\n", ans);
}

void kclist::listing(v_size deep, std::vector<v_size> & P, const std::vector<v_size> & C) {
    if(P.size() == k - 2) {
        for(ui i = 0; i < C.size(); i++) {
            for(ui j = i + 1; j < C.size(); j++) {
                if(cc(C[i], C[j])) ans++;
            }
        }

        return;
    }

    for(ui i = 0; i < C.size(); i++) {
        std::vector<v_size> newC;

        for(ui j = i + 1; j < C.size(); j++) {
            if(cc(C[i], C[j])) newC.push_back(C[j]);
        }

        P.push_back(C[i]);
        listing(deep + 1, P, newC);
        P.pop_back();
    }
}


#endif