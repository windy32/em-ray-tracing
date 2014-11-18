#include "Triangle.h"
#include <vector>
#include <math.h>

Triangle::Triangle()
{
    type = GeometryType::TRIANGLE;
}

Triangle::Triangle(const Point &a, const Point &b, const Point &c, const Vector &normal)
{
    this->a = a;
    this->b = b;
    this->c = c;
    this->normal = normal;
    type = GeometryType::TRIANGLE;
}

Triangle::Triangle(const Point &a, const Point &b, const Point &c)
{
    this->a = a;
    this->b = b;
    this->c = c;
    this->normal = Vector(a, b).cross(Vector(b, c)).norm();
    type = GeometryType::TRIANGLE;
}

static double det(double a11, double a12, double a13, 
                  double a21, double a22, double a23, 
                  double a31, double a32, double a33)
{
    return 
        a11 * a22 * a33 + 
        a12 * a23 * a31 + 
        a13 * a21 * a32 - 
        a13 * a22 * a31 - 
        a11 * a23 * a32 - 
        a12 * a21 * a33;
}

Point Triangle::getCenter() const
{
    return Point(
        (a.x + b.x + c.x) / 3,
        (a.y + b.y + c.y) / 3,
        (a.z + b.z + c.z) / 3);
}

void Triangle::getBoundingBox(Point &min, Point &max)
{
    double min_x = DBL_MAX, min_y = DBL_MAX, min_z = DBL_MAX;
    double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX; 

    if (a.x < min_x) min_x = a.x;
    if (b.x < min_x) min_x = b.x;
    if (c.x < min_x) min_x = c.x;

    if (a.x > max_x) max_x = a.x;
    if (b.x > max_x) max_x = b.x;
    if (c.x > max_x) max_x = c.x;

    if (a.y < min_y) min_y = a.y;
    if (b.y < min_y) min_y = b.y;
    if (c.y < min_y) min_y = c.y;

    if (a.y > max_y) max_y = a.y;
    if (b.y > max_y) max_y = b.y;
    if (c.y > max_y) max_y = c.y;

    if (a.z < min_z) min_z = a.z;
    if (b.z < min_z) min_z = b.z;
    if (c.z < min_z) min_z = c.z;

    if (a.z > max_z) max_z = a.z;
    if (b.z > max_z) max_z = b.z;
    if (c.z > max_z) max_z = c.z;

    min.x = min_x;
    min.y = min_y;
    min.z = min_z;

    max.x = max_x;
    max.y = max_y;
    max.z = max_z;
}

IntersectResult Triangle::intersect(Ray &ray)
{ 
    IntersectResult result(false);

    // A point P in triangle ABC:
    //  - P = alpha * A + beta * B + gamma * C
    //   (0 <= alpha <= 1, 0 <= beta <= 1, 0 <= gamma <= 1, and alpha + beta + gamma = 1)
    //  - P = (1 - beta - gamma) * A + beta * B + gamma * C
    //   (0 <= 1 - beta - gamma <= 1, 0 <= beta <= 1 and 0 <= gamma <= 1)
    //
    // Solve:
    //   O + t * DIR = (1 - beta - gamma) * A + beta * B + gamma * C
    //  ==> O + t * DIR = A + beta * (B - A) + gamma * (C - A)
    //  ==> (A - B) * beta + (A - C) * gamma + DIR * t = A - O
    //      [ x(A)-x(B)  x(A)-x(C)  x(DIR) ]   [ beta  ]   [ x(A) - x(O) ]
    //  ==> [ y(A)-y(B)  y(A)-y(C)  y(DIR) ] * [ gamma ] = [ y(A) - y(O) ]
    //      [ z(A)-z(B)  z(A)-z(C)  z(DIR) ]   [ t     ]   [ z(A) - z(O) ]
    //      [ m11  m12  m13 ]   [ beta  ]   [ b1 ]
    //  ==> [ m21  m22  m23 ] * [ gamma ] = [ b2 ]
    //      [ m31  m32  m33 ]   [ t     ]   [ b3 ]
    //  ==> M * [ beta  gamma  t ] ^ -1 = [ b1  b2  b3 ] ^ -1
    // 
    // According to the Cramer's Rule
    // http://en.wikipedia.org/wiki/Cramer%27s_rule
    //         | b1  m12  m13 |
    //  beta = | b2  m22  m23 | / | M |
    //         | b3  m32  m33 |
    //
    //          | m11  b1  m13 |
    //  gamma = | m21  b2  m23 | / | M |
    //          | m31  b3  m33 |
    //
    //      | m11  m12  b1 |
    //  t = | m21  m22  b2 | / | M |
    //      | m31  m32  b3 |
    double m11 = a.x - b.x;
    double m21 = a.y - b.y;
    double m31 = a.z - b.z;

    double m12 = a.x - c.x;
    double m22 = a.y - c.y;
    double m32 = a.z - c.z;

    double m13 = ray.direction.x;
    double m23 = ray.direction.y;
    double m33 = ray.direction.z;

    double b1 = a.x - ray.origin.x;
    double b2 = a.y - ray.origin.y;
    double b3 = a.z - ray.origin.z;

    double det_m = det(m11, m12, m13, m21, m22, m23, m31, m32, m33);
    if (fabs(det_m) < 1e-10)
    {
        return result;
    }

    double t = det(m11, m12, b1, m21, m22, b2, m31, m32, b3) / det_m;
    if (t < 0.0005f)
    {
        return result;
    }

    double beta = det(b1, m12, m13, b2, m22, m23, b3, m32, m33) / det_m;
    if (beta < -0.0001f || beta > 1.0001f) // avoid leaks
    {
        return result;
    }

    double gamma = det(m11, b1, m13, m21, b2, m23, m31, b3, m33) / det_m;
    if (gamma < -0.0001f || gamma > 1.0001f ||
        1 - beta - gamma < -0.0001f || 1 - beta - gamma > 1.0001f)
    {
        return result;
    }

    result.hit = true;
    result.geometry = this;
    result.distance = t;
    result.position = ray.getPoint(t);
    result.normal = normal;

    return result;
}

// Utils used by intersectWithGrid()
double getMin(const std::vector<Point> &points, Vector axis)
{
    double min = DBL_MAX; 
    for (unsigned int i = 0; i < points.size(); i++)
    {
        min = std::min(min, axis.dot(points[i]));
    }
    return min;
}

double getMax(const std::vector<Point> &points, Vector axis)
{
    double max = -DBL_MAX; 
    for (unsigned int i = 0; i < points.size(); i++)
    {
        max = std::max(max, axis.dot(points[i]));
    }
    return max;
}

bool intersectOnAxis(const std::vector<Point> &points1, const std::vector<Point> &points2, const Vector &axis)
{
    if (getMin(points1, axis) > getMax(points2, axis)) return false;
    if (getMax(points1, axis) < getMin(points2, axis)) return false;
    return true;     
}

// used for regular grid acceleration
bool Triangle::intersectWithGrid(const Grid &grid)
{
    std::vector<Point> gridPoints;
    gridPoints.push_back(grid.pos + Vector(0, 0, 0));
    gridPoints.push_back(grid.pos + Vector(0, 0, grid.size.z));
    gridPoints.push_back(grid.pos + Vector(0, grid.size.y, 0));
    gridPoints.push_back(grid.pos + Vector(0, grid.size.y, grid.size.z));
    gridPoints.push_back(grid.pos + Vector(grid.size.x, 0, 0));
    gridPoints.push_back(grid.pos + Vector(grid.size.x, 0, grid.size.z));
    gridPoints.push_back(grid.pos + Vector(grid.size.x, grid.size.y, 0));
    gridPoints.push_back(grid.pos + Vector(grid.size.x, grid.size.y, grid.size.z));

    std::vector<Point> trianglePoints;
    trianglePoints.push_back(a);
    trianglePoints.push_back(b);
    trianglePoints.push_back(c);

    // Test the x, y, and z axes
    if (!intersectOnAxis(gridPoints, trianglePoints, Vector(1, 0, 0))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, Vector(0, 1, 0))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, Vector(0, 0, 1))) return false;

    // Test the triangle normal
    if (!intersectOnAxis(gridPoints, trianglePoints, normal)) return false;

    // Test the 9 edge cross products
    Vector triangleEdge1(a, b);
    Vector triangleEdge2(b, c);
    Vector triangleEdge3(c, a);

    Vector boxEdge1 = Vector(1, 0, 0);
    Vector boxEdge2 = Vector(0, 1, 0);
    Vector boxEdge3 = Vector(0, 0, 1);

    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge1.cross(triangleEdge1))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge1.cross(triangleEdge2))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge1.cross(triangleEdge3))) return false;

    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge2.cross(triangleEdge1))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge2.cross(triangleEdge2))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge2.cross(triangleEdge3))) return false;

    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge3.cross(triangleEdge1))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge3.cross(triangleEdge2))) return false;
    if (!intersectOnAxis(gridPoints, trianglePoints, boxEdge3.cross(triangleEdge3))) return false;

    return true;
}
