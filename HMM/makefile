C++ = g++
CFLAGS1 = -Wall -g -c -std=c++11
CFLAGS2 = -g

all: hmm

hmm: hmm.o
	${C++} ${CFLAGS2} -o hmm hmm.o

hmm.o: hmm.cpp
	${C++} ${CFLAGS1} -o hmm.o hmm.cpp

clean:
	rm -f hmm hmm.o

