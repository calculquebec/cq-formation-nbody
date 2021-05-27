LIBS = -lm

CXX = pgc++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -fast -std=c++11 -Minfo=accel -acc -gpu=managed
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) #-DVERLET

nbody: nbody.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody nbody.cpp $(LIBS)
nbody_float: nbody_float.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody_float nbody_float.cpp $(LIBS)

clean:
	rm -f nbody nbody_float
