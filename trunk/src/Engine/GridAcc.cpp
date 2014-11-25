#include "GridAcc.h"
#include "Grid.h"
#include "Utils.h"
#include "Sphere.h"

std::vector<Geometry *> &GridAcc::get(int x, int y, int z)
{
    return data[(x * yLength + y) * zLength + z];
}

void GridAcc::getIndexInGrid(const Point &p, int &i, int &j, int&k)
{
    i = (int)((p.x - origin.x) / cellSizeX);
    j = (int)((p.y - origin.y) / cellSizeY);
    k = (int)((p.z - origin.z) / cellSizeZ);
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (k < 0) k = 0;
    if (i > xLength - 1) i = xLength - 1;
    if (j > yLength - 1) j = yLength - 1;
    if (k > zLength - 1) k = zLength - 1;
}

void GridAcc::init()
{
    Utils::PrintTime("Initialize Grid");

    // 1. Get the range of the triangles
    double min_x = DBL_MAX, min_y = DBL_MAX, min_z = DBL_MAX;
    double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX;

    for (unsigned int i = 0; i < scene->size(); i++)
    {
        Point min, max;
        (*scene)[i]->getBoundingBox(min, max);

        min_x = std::min(min_x, min.x);
        min_y = std::min(min_y, min.y);
        min_z = std::min(min_z, min.z);

        max_x = std::max(max_x, max.x);
        max_y = std::max(max_y, max.y);
        max_z = std::max(max_z, max.z);
    }

    double width = max_x - min_x;
    double height = max_y - min_y;
    double depth = max_z - min_z;
    double maxLength = std::max(std::max(width, height), depth);

    // Cut the longest dimension into 400 pieces
    double size = maxLength / 399;
    origin = Point(min_x - size / 2, min_y - size / 2, min_z - size / 2);
    cellSizeX = size;
    cellSizeY = size;
    cellSizeZ = size;
    xLength = (int)(width / cellSizeX + 1.5f);
    yLength = (int)(height / cellSizeY + 1.5f);
    zLength = (int)(depth / cellSizeZ + 1.5f);
    data.clear();
    data.resize(xLength * yLength * zLength);

    Utils::DbgPrint("Grid Size: %d x %d x %d\n", xLength, yLength, zLength);

    // For each triangle
    for (unsigned int m = 0; m < scene->size(); m++)
    {
        Point min, max;
        (*scene)[m]->getBoundingBox(min, max);

        int x_begin = (int)((min.x - origin.x) / cellSizeX);
        int y_begin = (int)((min.y - origin.y) / cellSizeY);
        int z_begin = (int)((min.z - origin.z) / cellSizeZ);

        int x_end = (int)((max.x - origin.x) / cellSizeX);
        int y_end = (int)((max.y - origin.y) / cellSizeY);
        int z_end = (int)((max.z - origin.z) / cellSizeZ);

        // Traverse the grid
        for (int i = x_begin; i <= x_end; i++)
        {
            for (int j = y_begin; j <= y_end; j++)
            {
                for (int k = z_begin; k <= z_end; k++)
                {
#if 0
                    Grid cell(grid.origin + Vector(i * size, j * size, k * size), Vector(size, size, size));
                    if (t->intersectWithGrid(cell))
                    {
                        grid.get(i, j, k).push_back(t);
                    }
#else
                    // Use a simple way to construct the grid, which is about 8 times faster.
                    // However, grid traversing is 20% slower
                    get(i, j, k).push_back((*scene)[m]);
#endif
                }
            }
        }
    }

#if 0 // Debug output
    for (int i = 0; i < grid.xLength; i++)
    {
        for (int j = 0; j < grid.yLength; j++)
        {
            Utils::DbgPrint("\n[%d, %d]", i, j);
            for (int k = 0; k < grid.zLength; k++)
            {
                Utils::DbgPrint(" %d", grid.get(i, j, k).size());
            }
        }
    }
#endif
}

IntersectResult GridAcc::intersect(Ray &ray, std::vector<RxIntersection> &rxPoints)
{
    Point near = origin;
    Point far = origin + Vector(
        cellSizeX * xLength, 
        cellSizeY * yLength, 
        cellSizeZ * zLength);

    // Current traversal state
    int cur_i, cur_j, cur_k; // the index in the grid
    double cur_d; // distance along the ray
    Point cur_p; // position

    // Is the origin of the ray outside of the grid?
    if (ray.origin.x < near.x || ray.origin.x > far.x ||
        ray.origin.y < near.y || ray.origin.y > far.y ||
        ray.origin.z < near.z || ray.origin.z > far.z)
    {
        double entryDistance, exitDistance;
        Grid sceneBox(origin, Vector(
            cellSizeX * xLength, 
            cellSizeY * yLength, 
            cellSizeZ * zLength));

        if (sceneBox.intersect(ray, entryDistance, exitDistance))
        {
            // Advance the ray to a grid boundary
            cur_d = entryDistance;
            cur_p = ray.getPoint(entryDistance);
            getIndexInGrid(cur_p, cur_i, cur_j, cur_k);
        }
        else
        {
            return IntersectResult(false);
        }
    }
    else // the origin of the ray is in the grid
    {
        cur_p = ray.origin;
        cur_d = 0;
        getIndexInGrid(ray.origin, cur_i, cur_j, cur_k);
    }

    std::map<int, std::pair<double, double> > rxIntersections;

    // Start traversing the grid
    while (true)
    {
        // See if the ray intersects with some triangle in the current cell
        std::vector<Geometry *> &list = get(cur_i, cur_j, cur_k);
        IntersectResult minResult(false);
        double minDistance = DBL_MAX;

        for (unsigned int i = 0; i < list.size(); i++)
        {
            IntersectResult result = list[i]->intersect(ray);
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
        
        if (minResult.hit)
        {
            std::map<int, std::pair<double, double> >::iterator it;
            for (it = rxIntersections.begin(); it != rxIntersections.end(); ++it)
            {
                if (it->second.first < minDistance)
                {
                    rxPoints.push_back(RxIntersection(it->first, it->second.first, it->second.second));
                }
            }

            return minResult; // There's a bug here, which should be fixed later
        }

        // Advance to the next cell with the 3D version of the DDA algorithm
        // http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
        Point p1 = origin + Vector(
            cur_i * cellSizeX, 
            cur_j * cellSizeY, 
            cur_k * cellSizeZ);
        Point p2 = origin + Vector(
            (cur_i + 1) * cellSizeX, 
            (cur_j + 1) * cellSizeY, 
            (cur_k + 1) * cellSizeZ);

        double dx = 0; // this is not delta_x, but is distance_on_x_axis
        double dy = 0; 
        double dz = 0;
        int di = 0;
        int dj = 0;
        int dk = 0;

        if (ray.direction.x > 0) // Plane x = p2.x
        {
            double cos_theta = ray.direction.x;
            dx = (p2.x - cur_p.x) / cos_theta;
            di = 1;
        }
        else // Plane x = p1.x
        {
            double cos_theta = -ray.direction.x;
            dx = (cur_p.x - p1.x) / cos_theta;
            di = -1;
        }

        if (ray.direction.y > 0) // Plane y = p2.y
        {
            double cos_theta = ray.direction.y;
            dy = (p2.y - cur_p.y) / cos_theta;
            dj = 1;
        }
        else // Plane y = p1.y
        {
            double cos_theta = -ray.direction.y;
            dy = (cur_p.y - p1.y) / cos_theta;
            dj = -1;
        }

        if (ray.direction.z > 0) // Plane z = p2.z
        {
            double cos_theta = ray.direction.z;
            dz = (p2.z - cur_p.z) / cos_theta;
            dk = 1;
        }
        else
        {
            double cos_theta = -ray.direction.z;
            dz = (cur_p.z - p1.z) / cos_theta;
            dk = -1;
        }

        // Advance
        if (dx < dy && dx < dz) // min = dx
        {
            cur_i += di;
            cur_d += dx;
            cur_p = ray.getPoint(cur_d);
        }
        else if (dy < dz) // min = dy
        {
            cur_j += dj;
            cur_d += dy;
            cur_p = ray.getPoint(cur_d);
        }
        else // min = dz
        {
            cur_k += dk;
            cur_d += dz;
            cur_p = ray.getPoint(cur_d);
        }

        // Leave the grid
        if (cur_i < 0 || cur_i > xLength - 1 ||
            cur_j < 0 || cur_j > yLength - 1 ||
            cur_k < 0 || cur_k > zLength - 1)
        {
            break;
        }
    }

    // Intersect with no triangles
    std::map<int, std::pair<double, double> >::iterator it;
    for (it = rxIntersections.begin(); it != rxIntersections.end(); ++it)
    {
        rxPoints.push_back(RxIntersection(it->first, it->second.first, it->second.second));
    }

    return IntersectResult(false);
}
