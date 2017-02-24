LIBS = -lboost_timer -lboost_system -lm

CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -march=native -funroll-loops -ffast-math
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) -DVERLET

nbody: nbody.cpp global.h
	$(CXX) $(CXX_FLAGS) -o nbody nbody.cpp $(LIBS)

omp: nbody_omp.cpp global.h
	$(CXX) $(CXX_FLAGS) -fopenmp -o nbody_omp nbody_omp.cpp $(LIBS)
    
clean:
	rm -f nbody











