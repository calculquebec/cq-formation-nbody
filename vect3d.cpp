#include <cmath>
#include <iomanip>

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


#define FNOP_VECT(op) Vect3d Vect3d::operator op(const Vect3d &vect) const \
    { return Vect3d(x op vect.x, \
                    y op vect.y, \
                    z op vect.z); }


#define FNOP_SCAL(op) Vect3d Vect3d::operator op(const double  scal) const \
    { return Vect3d(x op scal, \
                    y op scal, \
                    z op scal); }


#define FNOE_VECT(op) const Vect3d& Vect3d::operator op##=(const Vect3d &vect) \
    { x op##= vect.x; \
      y op##= vect.y; \
      z op##= vect.z; \
      return *this; }


#define FNOE_SCAL(op) const Vect3d& Vect3d::operator op##=(const double  scal) \
    { x op##= scal; \
      y op##= scal; \
      z op##= scal; \
      return *this; }


#define FNOP_LIST(op) FNOP_VECT(op) FNOP_SCAL(op) FNOE_VECT(op) FNOE_SCAL(op)


FNOP_LIST(+)
FNOP_LIST(-)
FNOP_LIST(*)
FNOP_LIST(/)


Vect3d floor(const Vect3d &v)
{
    return Vect3d(std::floor(v.x),
                  std::floor(v.y),
                  std::floor(v.z));
}


Vect3d operator*(double val, const Vect3d &v)
{
    return v * val;
}


std::ostream& operator<<(std::ostream &ostr, const Vect3d &v)
{
    for (int j = 0; j < 3; j++) {
        ostr << std::setw(10) << std::setprecision(4)
             << std::setiosflags(std::ios::fixed) << v.comp[j];
    }

    return ostr;
}


