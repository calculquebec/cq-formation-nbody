LIBS = -lm

CXX = pgc++
#CXX = icpc

DEBUG = -g -Wall

OPT = -O3 -fast -std=c++11 -Minfo=all,accel,intensity -ta:tesla:managed
#OPT = -O3 -xHost -ipo

CXX_FLAGS = $(OPT) #-DVERLET

nbody: nbody.cpp global.h Makefile
	$(CXX) $(CXX_FLAGS) -o nbody nbody.cpp $(LIBS)

clean:
	rm -f nbody











