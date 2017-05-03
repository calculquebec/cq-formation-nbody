#ifndef NBODY_H
#define NBODY_H


#include <vector>

#include "main.h"
#include "particule.h"


// Random number seed
#define DEFAULT_SEED 0
// Number of particles
#define DEFAULT_NP 100

// Time increment in ODE solver
#define DEFAULT_DT 0.05
// Number of time steps
#define DEFAULT_NT 10
// Frequency of writing position values to disk
#define DEFAULT_WRITE_FREQ 2
// Force softening
#define DEFAULT_EPSILON 0.00000000001

// Particle mass is drawn uniformly from the interval [low_mass, high_mass)
#define DEFAULT_LOW_MASS 1.0
#define DEFAULT_HIGH_MASS 5.0

// Fix the center of mass at 0
#define DEFAULT_CENTER_MASSES true
// Use a finite domain with toroidal boundary conditions?
#define DEFAULT_FINITE_DOMAIN false
// Condition for "bound state" of the particles
#define DEFAULT_BOUNDED_STATE true

// The dimensions of the finite domain
#define DEFAULT_L_MIN Vect3d(0.0, 0.0, 0.0)
#define DEFAULT_L_MAX Vect3d(100.0, 100.0, 50.0)


class NBody
{
public:
    enum Algo {Verlet, Runge_Kutta};

public:
    NBody(const Algo algo_ = Verlet);

    void configure(const Params &params);
    void integrate();

private:
    // Random number variables
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> VRG;

    // Verlet or Runge_Kutta
    Algo algo;

    // All particules
    std::vector<Particule> particules;

    // Time attributes: Delta T, N steps, write result freq, epsilon limit
    double dt;
    int NT;
    int write_freq;
    double epsilon;

    // Limits for the masses
    double low_mass;
    double high_mass;

    // Options
    bool center_masses;
    bool finite_domain;
    bool bounded_state;

    // Limits for the 3D space
    struct { Vect3d min; Vect3d max; } L;

};


#endif // NBODY_H


