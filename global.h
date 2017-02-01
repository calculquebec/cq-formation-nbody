// STL headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <string>
#include <limits>

// The Boost library headers...
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

// Number of particles
int NP = 100;
// Number of time steps
int NT = 10;
// Random number seed 
int seed = 0;
// Frequency of writing position values to disk
int write_freq = 2;
// Particle mass is drawn uniformly from the  
// interval [low_mass,high_mass)
double low_mass = 1.0;
double high_mass = 5.0;
// Time increment in ODE solver
double dt = 0.05;
// Force softening
double epsilon = 0.01;
// Use a finite domain with toroidal boundary 
// conditions?
bool finite_domain = false;
// The dimensions of the finite domain
double L[] = {0.0,100.0,0.0,100.0,0.0,50.0};

// Random number variables
boost::mt19937 BGT;
boost::uniform_real<> uniform;
boost::variate_generator<boost::mt19937&,boost::uniform_real<> >* VRG;


