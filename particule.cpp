#include "particule.h"


Particule::Particule()
{
}


void Particule::setMass(const double m)
{
    mass = m;
}


void Particule::setPos(const Vect3d &pos)
{
    position = pos;
}


void Particule::setVel(const Vect3d &vel)
{
    velocity = vel;
}


