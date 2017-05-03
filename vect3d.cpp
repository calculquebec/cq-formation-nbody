#include "vect3d.h"



Vect3d::Vect3d(): x(0.0), y(0.0), z(0.0)
{
}


Vect3d::Vect3d(double x_, double y_, double z_): x(x_), y(y_), z(z_)
{
}


Vect3d Vect3d::operator+(const Vect3d &v) const
{
    return Vect3d(x + v.x,
                  y + v.y,
                  z + v.z);
}


Vect3d Vect3d::operator/(double den) const
{
    return Vect3d(x / den,
                  y / den,
                  z / den);
}


