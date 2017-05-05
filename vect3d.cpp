#include "vect3d.h"


Vect3d::Vect3d()
{
}


Vect3d::Vect3d(double x_, double y_, double z_): x(x_), y(y_), z(z_)
{
}


double Vect3d::dotProd(const Vect3d &v) const
{
    return x * v.x +
           y * v.y +
           z * v.z;
}


Vect3d Vect3d::operator+(const Vect3d &v) const
{
    return Vect3d(x + v.x,
                  y + v.y,
                  z + v.z);
}


Vect3d Vect3d::operator-(const Vect3d &v) const
{
    return Vect3d(x - v.x,
                  y - v.y,
                  z - v.z);
}


Vect3d Vect3d::operator*(double val) const
{
    return Vect3d(x * val,
                  y * val,
                  z * val);
}


Vect3d Vect3d::operator/(double den) const
{
    return Vect3d(x / den,
                  y / den,
                  z / den);
}


const Vect3d& Vect3d::operator+=(const Vect3d &v)
{
    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}


const Vect3d& Vect3d::operator-=(const Vect3d &v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}


const Vect3d& Vect3d::operator*=(double val)
{
    x *= val;
    y *= val;
    z *= val;

    return *this;
}


const Vect3d& Vect3d::operator/=(double den)
{
    x /= den;
    y /= den;
    z /= den;

    return *this;
}


