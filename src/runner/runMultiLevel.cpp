#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../kclique/multilevel.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    
    string filePath = "data/skitter/";
    if(aC->exist("-f")) filePath = aC->get("-f");

    e_size N = 5000000;
    if(aC->exist("-N")) {
        string tmp = aC->get("-N");
        N = 0;
        for(long long unsigned i = 0; i < tmp.length(); i++) {
            N = N * 10 + (tmp[i]-'0');
        }
    }

    double alpha = 1.0;
    if(aC->exist("-a")) alpha = atof(aC->get("-a").c_str());

    v_size k = 10;
    if(aC->exist("-k")) k = atoi(aC->get("-k").c_str());

    Graph * g = new Graph();
    v_size n, maxK, tmp;
    double exCnt = 0.0;

    FILE * f = fopen((filePath+"s.txt").c_str(), "r");
    int err = fscanf(f, "%u", &n);
    if(err != 1) {
        printf("s.txt not exist\n");
        return  0;
    }
    if(~fscanf(f, "%u", &maxK)) {
        if(maxK >= k) {
            while(~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                if(tmp == k) break;
            }
        }
    }

    g->load(filePath+"edge.bin", filePath+"idx.bin", n);

    samplePlusExact * pt = new samplePlusExact(g);

    if(aC->exist("-cc")) 
        pt->runCC(k, deb, exCnt, alpha, N);
    else if(aC->exist("-cccpath"))
        pt->runCCPath(k, deb, exCnt, alpha, N);
    else if(aC->exist("-tripath"))
        pt->runTriPath(k, deb, exCnt, alpha, N);
    else
        pt->run(k, deb, exCnt, alpha, N);

    delete g;
    delete pt;
   

    

    return 0;
}