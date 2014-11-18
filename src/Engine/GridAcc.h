#ifndef GRID_ACC_H
#define GRID_ACC_H

#include "Accelerator.h"

class GridAcc : public Accelerator
{
private:
    Point origin;
    float cellSizeX;
    float cellSizeY;
    float cellSizeZ;
    int xLength;
    int yLength;
    int zLength;
    std::vector<std::vector<Triangle *>> data;

private:
    std::vector<Triangle *> &get(int x, int y, int z);
    void getIndexInGrid(const Point &p, int &i, int &j, int&k);

public:
    GridAcc(std::vector<Triangle *> *scene) : Accelerator(scene) {}
    virtual void init();
    virtual IntersectResult intersect(Ray &ray);
};

#endif
