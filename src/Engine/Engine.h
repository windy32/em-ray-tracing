#ifndef ENGINE_H
#define ENGINE_H

struct RtPoint
{
    double x;
    double y;
    double z;

    RtPoint() {}
    RtPoint(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct RtVector
{
    double x;
    double y;
    double z;

    RtVector() {}
    RtVector(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct RtTriangle
{
    RtPoint a;
    RtPoint b;
    RtPoint c;
    RtVector n;

    RtTriangle() {}
    RtTriangle(const RtPoint &a, const RtPoint &b, const RtPoint &c, const RtVector &n) : a(a), b(b), c(c), n(n) {}
};

enum RtPreprocessMethod
{
    Linear,
    Grid,
    KdTree
};

void Initialize();

void AddTriangle(const RtTriangle &triangle);
void AddTriangles(const RtTriangle *triangles, int n);
bool AddStlModel(const char *filename); // TODO: add unicode version

bool SetPreprocessMethod(RtPreprocessMethod method);
void SetTxPoint(const RtPoint &point, double power); // power in dBm
void SetRxPoints(const RtPoint *points, int n, double radius); // radius in meters

void SetParameters(
    double permittivity,
    double conductivity,
    int maxReflections,
    double raySpacing,   // unit: degree
    double frequency     // unit: MHz
    );

bool Simulate();
void GetRxPowers(double *powers, int n);

#endif
