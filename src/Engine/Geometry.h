#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Ray.h"
#include "IntersectResult.h"

enum GeometryType { TRIANGLE, SPHERE };

class Geometry
{
public:
    GeometryType type;

public:
    Geometry();
    virtual ~Geometry();
    virtual Point getCenter() const = 0;
    virtual void getBoundingBox(Point &min, Point &max) = 0;
    virtual IntersectResult intersect(Ray &ray) = 0;
};

#endif
