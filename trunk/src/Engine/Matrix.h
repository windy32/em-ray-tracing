#ifndef MATRIX_H
#define MATRIX_H

#include "Point.h"

class Matrix
{
public:
    double m11, m12, m13;
    double m21, m22, m23;
    double m31, m32, m33;

public:
    Matrix(
        double m11, double m12, double m13, 
        double m21, double m22, double m23,
        double m31, double m32, double m33);
    Point operator*(const Point &p) const;
    Vector operator*(const Vector &p) const;
};

#endif
