OBJS = main.o nbody.o particule.o vect3d.o
LIBS = -lm

CXX = g++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -funroll-loops -ffast-math -std=c++11
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) -DVERLET

nbody: $(OBJS)
	$(CXX) $(OPT) -o $@ $+ $(LIBS)

%.o: %.cpp %.h
	$(CXX) -c $(CXX_FLAGS) -o $@ $<

clean:
	rm -f nbody $(OBJS)
