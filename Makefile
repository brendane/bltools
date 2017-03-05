CXX = g++
CXXFLAGS = -I. --std=c++14 -Wall -O3
DEPS = SeqFileInWrapper.h

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

blhead: blhead.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o blhead blhead.o SeqFileInWrapper.o

blgrep: blgrep.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o blgrep blgrep.cpp SeqFileInWrapper.o

.PHONY: clean

clean:
	rm -f *.o blhead
