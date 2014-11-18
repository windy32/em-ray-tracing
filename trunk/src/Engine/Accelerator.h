#ifndef ACCELERATOR_H
#define ACCELERATOR_H

#include <vector>
#include "Triangle.h"

class Accelerator 
{
protected:
    std::vector<Triangle *> *scene;

public:
    Accelerator(std::vector<Triangle *> *scene) : scene(scene) {}
    virtual void init() = 0;
    virtual IntersectResult intersect(Ray &ray) = 0;
};

#endif
