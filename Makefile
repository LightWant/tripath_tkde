CXXFLAGS= -O3 -Wall -std=c++11 -fopenmp -march=native

all:bin/run bin/pivoter

bin/run:src/runner/run.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<

bin/lrun:src/runner/runLocal.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<

bin/makeTxt:src/runner/makeTxt.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/makeCSR:src/runner/makeCSR.cpp
	g++ $(CXXFLAGS) -o $@ $<
bin/sortData:src/runner/sortData.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/changeToD:src/runner/changeToD.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<
bin/tmp:src/runner/tmp.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<


bin/pivoter:src/runner/runPivoter.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/localColor:src/runner/runLocalColor.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/kclist:src/runner/kclist.cpp
	g++ $(CXXFLAGS) -o $@ $<

clean: 
	rm bin/*
