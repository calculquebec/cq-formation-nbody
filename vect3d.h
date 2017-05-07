#ifndef VECT3D_H
#define VECT3D_H


#include <ostream>


class Vect3d
{
public:
    Vect3d();
    Vect3d(double x_, double y_, double z_);

    double dotProd(const Vect3d &v) const;

    Vect3d operator+(const Vect3d &v) const;
    Vect3d operator-(const Vect3d &v) const;

    Vect3d operator*(double val) const;
    Vect3d operator/(double den) const;

    const Vect3d& operator+=(const Vect3d &v);
    const Vect3d& operator-=(const Vect3d &v);

    const Vect3d& operator*=(double val);
    const Vect3d& operator/=(double den);

public:
    union {
        struct { double x, y, z; };
        double comp[3];
    };
};


std::ostream& operator<<(std::ostream &ostr, const Vect3d &v);


#endif // VECT3D_H

