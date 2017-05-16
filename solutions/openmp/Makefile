LIBS = -lm

CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -march=native -funroll-loops -ffast-math -std=c++11 -fopenmp
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) #-DVERLET

nbody: nbody_omp.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody_omp nbody_omp.cpp $(LIBS)

clean:
	rm -f nbody











