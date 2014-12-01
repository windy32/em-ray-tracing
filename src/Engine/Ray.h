#ifndef RAY_H
#define RAY_H

#include "Vector.h"
#include <unordered_map>

#define COMPACT // compact path info to 4 bytes (unsigned int)

struct RayPath
{
#ifndef COMPACT
    int count; // number of reflection points
    int indexes[100]; // indexes of geometries (Geometry::index)
#endif
    int hash_code;

    static const unsigned int I = 17; // a small prime number
    static const unsigned int P = 486187739; // a big prime number

    RayPath();
    void addPoint(int index); // add a point to the reflection path
};

template <>
struct std::hash<RayPath>
{
    std::size_t operator()(const RayPath &key)
    {
        return key.hash_code;
    }
};

bool operator==(const RayPath &left, const RayPath &right);

struct Ray
{
    Point origin;
    Vector direction;

    // state & history information
    enum RayState { Start, FirstReflect, MoreReflect } state;
    double prev_mileage;
    Point  prev_point;

    // reflection path
    RayPath path;

    // ray tube
    double unit_surface_area;

    Ray(const Point &origin, const Vector &direction, double unitSurfaceArea) 
        : origin(origin), direction(direction), unit_surface_area(unitSurfaceArea),
          state(Start), prev_mileage(0), prev_point(origin)
    {
    }

    Point getPoint(double distance) const
    {
        return origin + direction * distance;
    }
};

#endif
