#ifndef ACCELERATOR_H
#define ACCELERATOR_H

#include <vector>
#include "Geometry.h"

class Accelerator 
{
protected:
    std::vector<Geometry *> *scene;

public:
    Accelerator(std::vector<Geometry *> *scene) : scene(scene) {}
    virtual void init() = 0;
    virtual IntersectResult intersect(Ray &ray) = 0;
};

#endif
