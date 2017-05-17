LIBS = -lm

CXX = nvcc
#CXX = icpc

DEBUG = -g -Wall

#OPT = -O3 -march=native -funroll-loops -ffast-math -std=c++11
OPT = -O3 -std=c++11
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) #-DVERLET

nbody: nbody.cu global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody nbody.cu $(LIBS)

clean:
	rm -f nbody











