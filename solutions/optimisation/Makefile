CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -march=native -funroll-loops -ffast-math -std=c++11
#OPT = -O3 -std=c++11 -xHost -ipo

CXX_FLAGS = -I../.. $(OPT) #-DVERLET

EXECUTABLES = nbody 1_nbody_transpose 2_nbody_permute_loops 3_nbody_sqrtf

all: $(EXECUTABLES)

nbody: nbody.cpp ../../global.h Makefile
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

1_nbody_transpose: 1_nbody_transpose.cpp ../../global.h Makefile
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

2_nbody_permute_loops: 2_nbody_permute_loops.cpp ../../global.h Makefile
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

3_nbody_sqrtf: 3_nbody_sqrtf.cpp ../../global.h Makefile
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(EXECUTABLES) *.o *.mol











