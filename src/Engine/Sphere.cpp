#include "Sphere.h"
#include <math.h>

Sphere::Sphere(const Point &center, double radius)
{
    this->center = center;
    this->radius = radius;
}

Point Sphere::getCenter() const
{
    return center;
}

void Sphere::getBoundingBox(Point &min, Point &max)
{
    min.x = center.x - radius;
    min.y = center.y - radius;
    min.z = center.z - radius;

    max.x = center.x + radius;
    max.y = center.y + radius;
    max.z = center.z + radius;
}

IntersectResult Sphere::intersect(Ray &ray)
{ 
    IntersectResult result(false);

    // Solve:
    //   | (o + t * dir) - c | = radius
    //  ==> ( co + t * dir ) ^ 2 - radius ^ 2 = 0
    //  ==> dir * dir * t ^ 2 + 2 * dir.co * t + (co.co - radius * radius) = 0
    //  ==> t ^ 2 + (2 * dir.co) * t + (co.co - radius * radius) = 0
    //  ==> t = - dir.co +- sqr((dir.co * dir.co) - (co.co - radius * radius))
    Vector co = Vector(center, ray.origin);
    double b = ray.direction.dot(co);
    double delta = b * b - (co.dot(co) - radius * radius);

    if (delta >= 0)
    {
        delta = sqrt(delta);
        if (-b + delta >= 0.0005f)
        {
            result.hit = true;
            result.geometry = this;
            result.distance = (-b - delta >= 0.0005f) ? -b - delta : -b + delta;
            result.position = ray.getPoint(result.distance);
            result.normal = Vector(center, result.position).norm();
        }
    }
    return result;
}

RxSphere::RxSphere(const Point &center, double radius, int index) : Sphere(center, radius)
{
    this->index = index;
}