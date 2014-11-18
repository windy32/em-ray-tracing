#include "GridAcc.h"
#include "KdTreeAcc.h"
#include "Triangle.h"
#include "Engine.h"

// Scene
std::vector<Geometry *> scene;

// Preprocessing
RtPreprocessMethod method;
GridAcc *accGrid = NULL;
KdTreeAcc *accKdTree = NULL;

// Tx point
Point txPoint;
double txPower;

// Rx points
std::vector<Point> rxPoints;
double rxRadius;

// Other parameters
struct RtParameter
{
    double permittivity;
    double conductivity;
    int maxReflections;
    double raySpacing;   // unit: degree
    double frequency;    // unit: MHz
} parameters;

void Initialize()
{
    scene.clear();
}

void AddTriangle(const RtTriangle &triangle)
{
    Point a = Point(triangle.a.x, triangle.a.y, triangle.a.z);
    Point b = Point(triangle.b.x, triangle.b.y, triangle.b.z);
    Point c = Point(triangle.c.x, triangle.c.y, triangle.c.z);
    Vector n = Vector(triangle.n.x, triangle.n.y, triangle.n.z);

    Triangle *t = new Triangle(a, b, c, n);
    scene.push_back(t);
}

void AddTriangles(const RtTriangle *triangles, int n)
{
    for (int i = 0; i < n; i++)
    {
        Point a = Point(triangles[i].a.x, triangles[i].a.y, triangles[i].a.z);
        Point b = Point(triangles[i].b.x, triangles[i].b.y, triangles[i].b.z);
        Point c = Point(triangles[i].c.x, triangles[i].c.y, triangles[i].c.z);
        Vector n = Vector(triangles[i].n.x, triangles[i].n.y, triangles[i].n.z);

        Triangle *t = new Triangle(a, b, c, n);
        scene.push_back(t);
    }
}

bool AddStlModel(const char *filename)
{
    // Open file
    FILE *fp = NULL;
    if (fopen_s(&fp, filename, "rb") != 0)
    {
        fprintf(stderr, "Error: Cannot open file \"%s\"\n", filename);
        return false;
    }

    // Read header and count
    char header[80];
    int count = -1;

    fread(header, 80, 1, fp);
    fread(&count, sizeof(int), 1, fp);

    // Read triangles
    for (int i = 0; i < count; i++)
    {
        float nx, ny, nz;
        float ax, ay, az;
        float bx, by, bz;
        float cx, cy, cz;
        short attribute;

        fread(&nx, sizeof(float), 1, fp);
        fread(&ny, sizeof(float), 1, fp);
        fread(&nz, sizeof(float), 1, fp);

        fread(&ax, sizeof(float), 1, fp);
        fread(&ay, sizeof(float), 1, fp);
        fread(&az, sizeof(float), 1, fp);

        fread(&bx, sizeof(float), 1, fp);
        fread(&by, sizeof(float), 1, fp);
        fread(&bz, sizeof(float), 1, fp);

        fread(&cx, sizeof(float), 1, fp);
        fread(&cy, sizeof(float), 1, fp);
        fread(&cz, sizeof(float), 1, fp);
        
        fread(&attribute, sizeof(short), 1, fp);

        Vector normal = Vector((double)nx, (double)ny, (double)nz);
        Point a = Point((double)ax, (double)ay, (double)az);
        Point b = Point((double)bx, (double)by, (double)bz);
        Point c = Point((double)cx, (double)cy, (double)cz);

        // Add to scene
        scene.push_back(new Triangle(a, b, c, normal));
    }

    return true;
}

bool Preprocess(RtPreprocessMethod m)
{
    method = m;
    if (method == Linear)
    {
        // Nothing to do here
    }
    else if (method == Grid)
    {
        accGrid = new GridAcc(&scene);
        // accGrid->init();
    }
    else if (method == KdTree)
    {
        accKdTree = new KdTreeAcc(&scene);
        // accKdTree->init();
    }
    else
    {
        fprintf(stderr, "Error: Unknown preprocess method\n");
        return false;
    }

    return true;
}

void SetTxPoint(const RtPoint &point, double power)
{
    txPoint = Point(point.x, point.y, point.z);
    txPower = power;
}

void SetRxPoints(const RtPoint *points, int n, double radius)
{
    rxPoints.clear();
    for (int i = 0; i < n; i++)
    {
        rxPoints.push_back(Point(points[i].x, points[i].y, points[i].z));
    }
    rxRadius = radius;
}

void SetParameters(
    double permittivity, double conductivity,  int maxReflections,
    double raySpacing, double frequency)
{
    parameters.permittivity = permittivity;
    parameters.conductivity = conductivity;
    parameters.maxReflections = maxReflections;
    parameters.raySpacing = raySpacing;
    parameters.frequency = frequency;
}

bool Simulate()
{
    // TODO: print warning messages
    //       when part of the parameters have not been specified
    return true;
}

bool GetRxPowers(double *powers, int n)
{
    return true;
}
