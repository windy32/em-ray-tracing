#ifndef COMPLEX_H
#define COMPLEX_H

#include "Vector.h"

class ComplexVector;

class ComplexNumber
{
public:       // a + bi
    double a; // real part
    double b; // imaginary part

public:
    ComplexNumber(double a, double b);

    ComplexNumber operator+(const ComplexNumber &v) const;
    ComplexNumber operator-(const ComplexNumber &v) const;
    ComplexNumber operator+(double v) const;
    ComplexNumber operator-(double v) const;

    ComplexNumber operator*(const ComplexNumber &v) const;
    ComplexNumber operator*(double v) const;
    ComplexVector operator*(const Vector &v) const;

    ComplexNumber operator/(const ComplexNumber &v) const;

    ComplexNumber Sqrt() const;

    // Calculate A * e ^ (i * phi)
    //         = A * cos(phi) * A * i * sin(phi)
    static ComplexNumber Euler(double A, double phi);
};

class ComplexVector
{
public:
    ComplexNumber x;
    ComplexNumber y;
    ComplexNumber z;

public:
    ComplexVector(ComplexNumber x, ComplexNumber y, ComplexNumber z);

    ComplexVector operator+(const ComplexVector &v);
    ComplexVector operator*(double v);
};

#endif
