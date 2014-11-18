#ifndef GRID_ACC_H
#define GRID_ACC_H

#include "Accelerator.h"

class GridAcc : public Accelerator
{
private:
    Point origin;
    double cellSizeX;
    double cellSizeY;
    double cellSizeZ;
    int xLength;
    int yLength;
    int zLength;
    std::vector<std::vector<Geometry *>> data;

private:
    std::vector<Geometry *> &get(int x, int y, int z);
    void getIndexInGrid(const Point &p, int &i, int &j, int&k);

public:
    GridAcc(std::vector<Geometry *> *scene) : Accelerator(scene) {}
    virtual void init();
    virtual IntersectResult intersect(Ray &ray);
};

#endif
