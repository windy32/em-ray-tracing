#ifndef VECTOR_H
#define VECTOR_H

#include "Point.h"

#include <math.h>

const double PI = 3.14159265359;

class Point;

class Vector
{
public:
    double x;
    double y;
    double z;

public:
    Vector(double x = 0, double y = 0, double z = 0);
    Vector(const Point &start, const Point &end);

    double length() const;
    double sqrLength() const;
    double & operator[](int index);
    const double & operator[](int index) const;

    Vector operator+(const Vector &b) const;
    Vector operator-(const Vector &b) const;
    Vector operator*(double b) const;
    Vector mult(const Vector &b) const;
    Vector& norm();
    double dot(const Vector &b) const;
    double dot(const Point &p) const;
    Vector cross(const Vector &b);

    double angleTo(const Vector &b) const;
};

#endif
