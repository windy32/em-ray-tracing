#ifndef RAY_H
#define RAY_H

#include "Vector.h"

struct Ray
{
    Point origin;
    Vector direction;

    // state & history information
    enum RayState { Start, FirstReflect, MoreReflect } state;
    double prev_mileage;
    Point  prev_point;

    Ray(const Point &origin, const Vector &direction) 
        : origin(origin), direction(direction), state(Start), prev_mileage(0), prev_point(origin)
    {
    }

    Point getPoint(double distance) const
    {
        return origin + direction * distance;
    }
};

#endif
