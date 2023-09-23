#include "../tools/getArgs.hpp"
#include "../graph/graph.hpp"
#include "../kclique/kclist.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;
/**
 * nvcc src\runPivoter.cpp -rdc=true --gpu-architecture=sm_70 -I "D:\cuda\NIVIDIA  GPU computing Toolklt_v10.1\CUDA\v10.1\include" -I "D:\cuda\NVIDIA Corporation_v10.1\common\inc" -o bin\runPivoter src\kClique\pivoterGPU.cu  -Xcompiler "/wd 4819"
 * **/
int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);

    string filePath = "data/skitter/";
    if(aC->exist("-f")) filePath = aC->get("-f");

    v_size k = 5;
    if(aC->exist("-k")) k = atoi(aC->get("-k").c_str());

    std::cout << "file:" << filePath << std::endl;
    std::cout << "k:" << k << std::endl;

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
    printf("degeneracy %u\n", g->degeneracy);

    kclist * pt = new kclist(g, k);
    
    auto s1 = std::chrono::steady_clock::now();
    // pt->run(exCnt);
    pt->runLinearset(exCnt);
    auto s2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
    std::cout << "time:" << duration.count() << "ms" << std::endl;

    delete pt;
    delete g;

    return 0;
}