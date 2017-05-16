CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -march=native -funroll-loops -ffast-math -std=c++11
#OPT = -O3 -std=c++11 -xHost -ipo

CXX_FLAGS = -I../.. $(OPT) #-DVERLET

EXECUTABLES = nbody_opt_omp

all: $(EXECUTABLES)

nbody: nbody.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

nbody_opt_omp: nbody_opt_omp.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -fopenmp -o $@ $< $(LIBS)

clean:
	rm -f $(EXECUTABLES) *.o *.mol











