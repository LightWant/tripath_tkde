#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../localClique/localCounter.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    
    v_size deb = 0;
    if(aC->exist("-deb")) deb = atoi(aC->get("-deb").c_str());
    
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
printf("degeneracy %u\n", g->degeneracy);fflush(stdout);

    localCounter * pt = new localCounter(g);

    auto ck = [&](vector<v_size> & ans, vector<v_size> & nodes) {
        std::vector<v_size> real(g->vCnt);
        std::ifstream fin(filePath+aC->get("-k")+"local.txt");
        for(ui i = 0; i < g->vCnt; i++) fin >> real[i];
        fin.close();

        std::vector<v_size> errs(10);
        v_size largerThan10 = 0, largerThan10And100 = 0;
        
        for(ui u : nodes) {
            ui i = u;
            // double err = std::abs((double)real[i] - (double)ans[i]) / real[i];

            if(real[i] == 0) {
                if(ans[i] > 10) {
                    largerThan10++;
                    // std::cout <<  i << " " << 
                    //     real[i] << ' ' << ans[i] << std::endl; 
                }
                else  {
                    errs[1]++;
                }
                continue;
            }
            if(ans[i] == 0) {
                if(real[i] > 10) {
                    largerThan10++;
                    // std::cout <<  i << " " << 
                    //     real[i] << ' ' << ans[i] << std::endl; 
                } 
                else {
                    errs[1]++;
                }
                continue;
            }
            double err = std::max((double)real[i]/ans[i], (double)ans[i]/real[i]);
            
            if(std::round(err) < 10) {
                errs[std::round(err)]++;
            }
            else {
                largerThan10++;
                // if(real[i] > 100) largerThan10And100++;
                // std::cout <<  i << " " << err << ' ' << 
                //     real[i] << ' ' << ans[i] << std::endl; 
            }
            // if(err > 2)
            //     std::cout <<  i << " " << err << ' ' << 
            //     real[i] << ' ' << ans[i] << std::endl; 
        }

        for(ui i = 1; i < 10; i++) {
            std::cout << "err_"<<i<<":" << errs[i] << std::endl;
        }
        std::cout << "errlarger10:" << largerThan10 << std::endl;
        // std::cout << "errlarger10And100:" << largerThan10And100 << std::endl;
    };

    if(aC->exist("-cc")) 
        pt->runCC(k, deb, exCnt, alpha, N);
    else if(aC->exist("-cccpath"))
        pt->runCCPath(k, deb, exCnt, alpha, N);
    else if(aC->exist("-tripath")) {
        std::vector<v_size> ans = pt->runTriPath(k, deb, exCnt, alpha, N);
std::vector<v_size> nodes = pt->getNodes();
printf("node size %u\n", nodes.size());fflush(stdout);
        if(aC->exist("-check")) ck(ans, nodes);
    }
    else if(aC->exist("-allccpath")) {
        std::vector<v_size> ans = pt->runAllCCPath(k, deb, N);
std::vector<v_size> nodes = pt->getNodes();
        if(aC->exist("-check")) ck(ans, nodes);
    }
    else if(aC->exist("-all"))
        pt->countExactAll(k, filePath+aC->get("-k")+"local.txt");

    delete g;
    delete pt;
    

    return 0;
}