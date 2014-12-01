#include "../Engine/Engine.h"
#include <stdio.h>
#include <math.h>

#define N 300

void freespace()
{
    // 1. init
    Initialize();

    // 2. do not add any triangle

    // 3. set preprocessing method
    SetPreprocessMethod(RtPreprocessMethod::KdTree);

    // 4. set tx point
    SetTxPoint(RtPoint(0, 0, 0), 20);

    // 4. set rx points
    RtPoint rxPoints[N];
    for (int i = 0; i < N; i++)
    {
        rxPoints[i].x = i + 1;
        rxPoints[i].y = 0;
        rxPoints[i].z = 0;
    }
    SetRxPoints(rxPoints, N, 0.6545);

    // 5. set parameters
    SetParameters(
        7.0,     // relative permittivity of concrete
        0.0015,  // electrical conductivity of concrete
        5,       // max reflectinos
        0.25,    // ray spacing (degrees)
        2437.0); // frequency (MHz)

    // 6. simulate
    Simulate();

    // 7. get results
    double rxPowers[N];
    GetRxPowers(rxPowers, N);

    // 8. output
    for (int i = 0; i < N; i++)
    {
        printf("%.2lf\n", rxPowers[i]);
    }
}

void twoRayGround()
{
    // 1. init
    Initialize();

    // 2. do not add any triangle
    AddTriangle(RtTriangle(
        RtPoint(0, -N, 0),
        RtPoint(0, N, 0),
        RtPoint(N, 0, 0),
        RtVector(0, 0, 1)));

    // 3. set preprocessing method
    SetPreprocessMethod(RtPreprocessMethod::KdTree);

    // 4. set tx point
    SetTxPoint(RtPoint(0, 0, 2), 20);

    // 4. set rx points
    RtPoint rxPoints[N];
    for (int i = 0; i < N; i++)
    {
        rxPoints[i].x = i + 1;
        rxPoints[i].y = 0;
        rxPoints[i].z = 2;
    }
    SetRxPoints(rxPoints, N, 0.6545);

    // 5. set parameters
    SetParameters(
        7.0,     // relative permittivity of concrete
        0.0015,  // electrical conductivity of concrete
        1,       // max reflectinos
        0.25,    // ray spacing (degrees)
        2437.0); // frequency (MHz)

    // 6. simulate
    Simulate();

    // 7. get results
    double rxPowers[N];
    GetRxPowers(rxPowers, N);

    // 8. output
    for (int i = 0; i < N; i++)
    {
        printf("%.2lf\n", rxPowers[i]);
    }
}

void tunnel()
{
    // 1. init
    Initialize();

    // 2. do not add any triangle
    if (!AddStlModel("tunnel-60-90.stl"))
        return;

    // 3. set preprocessing method
    SetPreprocessMethod(RtPreprocessMethod::KdTree);

    // 4. set tx point
    SetTxPoint(RtPoint(0, 0, 4), 20);

    // 4. set rx points
    RtPoint rxPoints[300];
    for (int i = 0; i < 300; i++)
    {
        double theta = (i + 1) * 3.14159265358979324 * 0.5 / 300.0;
        rxPoints[i].x = 300 * (1.0 - cos(theta));
        rxPoints[i].y = 300 * sin(theta);
        rxPoints[i].z = 4.0;
    }
    SetRxPoints(rxPoints, 300, 0.5);

    // 5. set parameters
    SetParameters(
        7.0,     // relative permittivity of concrete
        0.0015,  // electrical conductivity of concrete
        8,       // max reflectinos
        0.25,    // ray spacing (degrees)
        2437.0); // frequency (MHz)

    // 6. simulate
    Simulate();

    // 7. get results
    double rxPowers[300];
    GetRxPowers(rxPowers, 300);

    // 8. output
    for (int i = 0; i < 300; i++)
    {
        printf("%.2lf\n", rxPowers[i]);
    }
}

int main(int argc, char *argv[])
{
    twoRayGround();
    return 0;
}
