#ifndef ACCELERATOR_H
#define ACCELERATOR_H

#include <vector>
#include "Geometry.h"

struct RxSphereInfo // Used by accelerators
{
    double distance;
    double offset;
    double radius;
};

class Accelerator 
{
protected:
    std::vector<Geometry *> *scene;

public:
    Accelerator(std::vector<Geometry *> *scene) : scene(scene) {}
    virtual void init() = 0;
    virtual IntersectResult intersect(Ray &ray, std::vector<RxIntersection> &rxPoints) = 0;
};

#endif
