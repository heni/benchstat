run-benchmarks: run-benchmarks.cpp
	g++ -O2 -std=c++0x $< -o $@
