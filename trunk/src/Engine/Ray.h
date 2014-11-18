#ifndef RAY_H
#define RAY_H

#include "Vector.h"

struct Ray
{
    Point origin;
    Vector direction;

    Ray(const Point &origin, const Vector &direction) : origin(origin), direction(direction)
    {
    }

    Point getPoint(float distance) const
    {
        return origin + direction * distance;
    }
};

#endif
