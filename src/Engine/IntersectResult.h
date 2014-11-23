#ifndef INTERSECT_RESULT_H
#define INTERSECT_RESULT_H

#include "Vector.h"
#include "Geometry.h"
#include <map>

class Geometry;

//        rx sphere 1   rx sphere 2    /
//           +----+      +----+       /
// --------->|    |----->|    |----->/ intersection point
//  ray      +----+      +----+     /
//                                 /
struct RxIntersection
{
    double distance;
    int index;

    RxIntersection(int index, double distance)
    {
        this->index = index;
        this->distance = distance;
    }
};

struct IntersectResult
{
    bool      hit;
    Geometry* geometry;
    double    distance;
    Point     position;

    // The normal vector that points to the outside of the object
    Vector    normal;

    IntersectResult()
    {
    }

    IntersectResult(bool hit)
    {
        this->hit = hit; 
    }
};

#endif
