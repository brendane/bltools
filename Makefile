CXX = g++
CXXFLAGS = -I. --std=c++14 -Wall -O3 -Ideps/tclap/include

blgrep: blgrep.cpp
	$(CXX) $(CXXFLAGS) -o blgrep blgrep.cpp
