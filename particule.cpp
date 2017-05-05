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


const double& Particule::m() const
{
    return mass;
}


double& Particule::m()
{
    return mass;
}


const Vect3d& Particule::p() const
{
    return position;
}


Vect3d& Particule::p()
{
    return position;
}


const Vect3d& Particule::v() const
{
    return velocity;
}


Vect3d& Particule::v()
{
    return velocity;
}


