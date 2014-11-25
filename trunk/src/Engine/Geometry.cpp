#include "Geometry.h"

#define NULL 0

Geometry::Geometry()
{
    index = ++count;
}

Geometry::~Geometry()
{
}

int Geometry::count = 0;