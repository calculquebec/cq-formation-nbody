#ifndef PARTICULE_H
#define PARTICULE_H


#include "vect3d.h"


class Particule
{
public:
    Particule();

    void setMass(const double m);
    void setPos(const Vect3d &pos);
    void setVel(const Vect3d &vel);

private:
    double mass;
    Vect3d position;
    Vect3d velocity;

};


#endif	// PARTICULE_H


