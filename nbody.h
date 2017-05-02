#ifndef NBODY_H
#define NBODY_H


#include <vector>

#include "main.h"
#include "particule.h"


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
    bool finite_domain;
    bool center_masses;
    bool bounded_state;

    // Limits for the 3D space
    struct { Vect3d min; Vect3d max; } L;

};


#endif // NBODY_H


