#include "Ray.h"

RayPath::RayPath()
{
#ifndef COMPACT
    count = 0;
    memset(indexes, 0, 100 * sizeof(int));
#endif
    hash_code = I;
}

void RayPath::addPoint(int index)
{
#ifndef COMPACT
    indexes[count++] = index;
#endif
    hash_code = hash_code * P + index;
}

bool operator==(const RayPath &left, const RayPath &right)
{
#ifndef COMPACT
    if (left.count != right.count)
        return false;

    for (int i = 0; i < left.count; i++)
    {
        if (left.indexes[i] != right.indexes[i])
            return false;
    }
    return true;
#else
    return (left.hash_code == right.hash_code);
#endif
}
