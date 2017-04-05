// STL headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <random>
#include <limits>

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
double epsilon = 0.00000000001;
// Use a finite domain with toroidal boundary 
// conditions?
bool finite_domain = false;
// The dimensions of the finite domain
double L[] = {0.0,100.0,0.0,100.0,0.0,50.0};
// Condition for "bound state" of the particles
bool bounded_state = true;
// Fix the center of mass at 0
bool center_masses = true;

// Random number variables
std::random_device rd;  
std::mt19937 gen(rd());
std::uniform_real_distribution<> VRG(0.0,1.0);


