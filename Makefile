CXX = g++
CXXFLAGS = -I. --std=c++14 -Wall -O3 -fPIC
DEPS = SeqFileInWrapper.h

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

blwc: blwc.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o blwc blwc.o SeqFileInWrapper.o

blhead: blhead.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o blhead blhead.o SeqFileInWrapper.o

bltail: bltail.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o bltail bltail.o SeqFileInWrapper.o

blgrep: blgrep.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o blgrep blgrep.cpp SeqFileInWrapper.o

bljoin: bljoin.o SeqFileInWrapper.o
	$(CXX) $(CXXFLAGS) -o bljoin bljoin.cpp SeqFileInWrapper.o

.PHONY: clean

clean:
	rm -f *.o blhead bltail blwc blgrep bljoin
