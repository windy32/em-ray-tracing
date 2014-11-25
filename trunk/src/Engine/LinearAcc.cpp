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
    //    - Value 1 (first): distance - distance from the origin of the ray to the intersection point
    //    - Value 2 (second): offset - distance from the intersection point to the center of rx sphere
    // Note: A ray may intersect with rx spheres even when result.hit = false
    std::map<int, std::pair<double, double> > rxIntersections;

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
                rxIntersections[s->index].first = result.distance;
                rxIntersections[s->index].second = Vector(result.position, s->center).length();
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

    std::map<int, std::pair<double, double> >::iterator it;
    for (it = rxIntersections.begin(); it != rxIntersections.end(); ++it)
    {
        if (it->second.first < minDistance)
        {
            rxPoints.push_back(RxIntersection(it->first, it->second.first, it->second.second));
        }
    }

    return minResult;
}
