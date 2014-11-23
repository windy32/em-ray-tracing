#include "../Engine/Engine.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    // 1. init
    Initialize();

    // 2. do not add any triangle

    // 3. set preprocessing method
    SetPreprocessMethod(RtPreprocessMethod::Linear);

    // 4. set tx point
    SetTxPoint(RtPoint(0, 0, 0), 20);

    // 4. set rx points
    RtPoint rxPoints[300];
    for (int i = 0; i < 300; i++)
    {
        rxPoints[i].x = 0;
        rxPoints[i].y = i + 1;
        rxPoints[i].z = 0;
    }
    SetRxPoints(rxPoints, 300, 0.6545);

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
    double rxPowers[300];
    GetRxPowers(rxPowers, 300);

    // 8. output
    for (int i = 0; i < 300; i++)
    {
        printf("%.2lf\n", rxPowers[i]);
    }

    return 0;
}
