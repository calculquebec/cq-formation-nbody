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
    Algo algo;

    std::vector<Particule> particules;

};


#endif // NBODY_H


