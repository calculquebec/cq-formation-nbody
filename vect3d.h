#ifndef VECT3D_H
#define VECT3D_H


#include <ostream>


class Vect3d
{
public:
    Vect3d();
    Vect3d(double x_, double y_, double z_);

    double dotProd(const Vect3d &v) const;

    #define OP_VECT(op)       Vect3d  operator op   (const Vect3d &vect) const;
    #define OP_SCAL(op)       Vect3d  operator op   (const double  scal) const;
    #define OE_VECT(op) const Vect3d& operator op##=(const Vect3d &vect)      ;
    #define OE_SCAL(op) const Vect3d& operator op##=(const double  scal)      ;

    #define OP_LIST(op) OP_VECT(op) OP_SCAL(op) OE_VECT(op) OE_SCAL(op)

    OP_LIST(+)
    OP_LIST(-)
    OP_LIST(*)
    OP_LIST(/)

public:
    union {
        struct { double x, y, z; };
        double comp[3];
    };
};


Vect3d floor(const Vect3d &v);
Vect3d operator*(double val, const Vect3d &v);
std::ostream& operator<<(std::ostream &ostr, const Vect3d &v);


#endif // VECT3D_H

