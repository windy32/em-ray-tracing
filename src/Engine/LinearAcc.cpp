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
    //    - Value: distance
    // Note: A ray may intersect with rx spheres even when result.hit = false
    std::map<int, double> rxIntersections;

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
                rxIntersections[s->index] = result.distance;
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

    std::map<int, double>::iterator it;
    for (it = rxIntersections.begin(); it != rxIntersections.end(); ++it)
    {
        if (it->second < minDistance)
        {
            rxPoints.push_back(RxIntersection(it->first, it->second));
        }
    }

    return minResult;

}
