#ifndef POINT_H
#define POINT_H

#include "Vector.h"

class Vector;

class Point
{
public:
    double x;
    double y;
    double z;

public:
    Point(double x = 0, double y = 0, double z = 0);

    Point operator+(const Vector &v) const;
    double & operator[](int index);
    const double & operator[](int index) const;
};

#endif
