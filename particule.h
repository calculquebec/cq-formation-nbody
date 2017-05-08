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
    void setAcc(const Vect3d &acc);

    const double& m() const;
    double& m();

    const Vect3d& p() const;
    Vect3d& p();

    const Vect3d& v() const;
    Vect3d& v();

    const Vect3d& a() const;
    Vect3d& a();

    double potentialEnergy(const Particule &part, const double epsilon = 0.0) const;
    double kineticEnergy() const;

private:
    double mass;
    Vect3d position;
    Vect3d velocity;
    Vect3d acceleration;

};


#endif	// PARTICULE_H


