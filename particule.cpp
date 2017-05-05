#include <cmath>

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


double Particule::potentialEnergy(const Particule &part, const double epsilon) const
{
    Vect3d delta = position - part.position;
    return mass * part.mass / std::sqrt(epsilon + delta.dotProd(delta));
}


double Particule::kineticEnergy() const
{
    return 0.5 * mass * velocity.dotProd(velocity);
}


