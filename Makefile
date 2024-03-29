LIBS = -lm

CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -funroll-loops -ffast-math -std=c++11
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT)  # -DVERLET

nbody: nbody.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody nbody.cpp $(LIBS)

clean:
	rm -f nbody
