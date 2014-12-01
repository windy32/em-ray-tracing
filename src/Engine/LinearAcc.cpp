#include "LinearAcc.h"
#include "Sphere.h"

void LinearAcc::init()
{
    // nothing to do here
}

IntersectResult LinearAcc::intersect(Ray &ray, std::vector<RxIntersection> &rxPoints)
{
    // Intersection points with rx spheres
    //    - Key: rx sphere index
    //    - Value: rx sphere info (distance, offset, radius)
    //       - distance: distance from the origin of the ray to the intersection point
    //       - offset: distance from the intersection point to the center of rx sphere
    // Note: A ray may intersect with rx spheres even when result.hit = false
    std::map<int, RxSphereInfo> rxIntersections;

    double minDistance = DBL_MAX;
    IntersectResult minResult(false);

    for (unsigned int i = 0; i < scene->size(); i++)
    {
        IntersectResult result = (*scene)[i]->intersect(ray);
        if (result.hit)
        {
            if (result.geometry->type == SPHERE && // rx sphere
                result.distance < minDistance)
            {
                RxSphere *s = (RxSphere *)result.geometry;
                rxIntersections[s->index].distance = result.distance;
                rxIntersections[s->index].offset = Vector(result.position, s->center).length();
                rxIntersections[s->index].radius = s->radius;
            }
            else // triangle
            {
                if (result.distance < minDistance) 
                {
                    minDistance = result.distance;
                    minResult = result;
                }
            }
        }
    }

    std::map<int, RxSphereInfo>::iterator it;
    for (it = rxIntersections.begin(); it != rxIntersections.end(); ++it)
    {
        if (it->second.distance < minDistance)
        {
            rxPoints.push_back(
                RxIntersection(it->first, it->second.distance, it->second.offset, it->second.radius));
        }
    }

    return minResult;
}
