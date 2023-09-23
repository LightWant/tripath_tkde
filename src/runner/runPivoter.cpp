#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

// #include "../kClique/ccpath.hpp"
#include "../kclique/pivoterMsk.hpp"
#include "../kclique/pivoterX.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;
#include <cassert>
/**
 * nvcc src\runPivoter.cpp -rdc=true --gpu-architecture=sm_70 -I "D:\cuda\NIVIDIA  GPU computing Toolklt_v10.1\CUDA\v10.1\include" -I "D:\cuda\NVIDIA Corporation_v10.1\common\inc" -o bin\runPivoter src\kClique\pivoterGPU.cu  -Xcompiler "/wd 4819"
 * **/
int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    Graph * g = new Graph();
    g->load(aC->get("-edge"), aC->get("-idx"), atoi(aC->get("-v").c_str()));

    printf("degeneracy %u\n", g->degeneracy);fflush(stdout);

    PivoterMsk * pt = new PivoterMsk(g);
    // PivoterX * pt = new PivoterX(g);
    v_size deb = atoi(aC->get("-debug").c_str());
    // if(aC->get("-debug").c_str() != "") deb = ;
    int k = 5;
    if(aC->exist("-k")) k = atoi(aC->get("-k").c_str());
    pt->runV(k, deb);

    delete pt;
    delete g;
    delete aC;

    return 0;
}