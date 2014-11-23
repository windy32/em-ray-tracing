#ifndef LINEAR_ACC_H
#define LINEAR_ACC_H

#include "Accelerator.h"

class LinearAcc : public Accelerator
{
public:
    LinearAcc(std::vector<Geometry *> *scene) : Accelerator(scene) {}
    virtual void init();
    virtual IntersectResult intersect(Ray &ray, std::vector<RxIntersection> &rxPoints);
};

#endif
