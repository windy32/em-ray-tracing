#include "Complex.h"
#include <math.h>

ComplexNumber::ComplexNumber(double a, double b)
{
    this->a = a;
    this->b = b;
}

ComplexNumber ComplexNumber::operator+(const ComplexNumber &v) const
{
    return ComplexNumber(a + v.a, b + v.b);
}

ComplexNumber ComplexNumber::operator-(const ComplexNumber &v) const
{
    return ComplexNumber(a - v.a, b - v.b);
}

ComplexNumber ComplexNumber::operator+(double v) const
{
    return ComplexNumber(a + v, b);
}

ComplexNumber ComplexNumber::operator-(double v) const
{
    return ComplexNumber(a - v, b);
}

ComplexNumber ComplexNumber::operator*(const ComplexNumber &v) const
{
    // (a + bi)(c + di) = (ac - bd) + (bc + ad)i
    double c = v.a;
    double d = v.b;
    return ComplexNumber(a * c - b * d, b * c + a * d);
}

ComplexNumber ComplexNumber::operator*(double v) const
{
    return ComplexNumber(a * v, b * v);
}

ComplexVector ComplexNumber::operator*(const Vector &v) const
{
    return ComplexVector((*this) * v.x, (*this) * v.y, (*this) * v.z);
}

ComplexNumber ComplexNumber::operator/(const ComplexNumber &v) const
{
    // http://en.wikipedia.org/wiki/Complex_number
    //
    // a + bi   (a + bi)*(c - di)      ac + bd        bc - ad
    // ------ = ----------------- = (----------) + (----------)i
    // c + di   (c + di)*(c - di)     c^2 + d^2       c^2 + d^2
    double c = v.a;
    double d = v.b;
    double r = sqrt(c * c + d * d);
    return ComplexNumber((a * c + b * d) / r, (b * c - a * d) / r);
}

ComplexNumber ComplexNumber::Sqrt() const
{
    // The square roots of a + bi (b != 0) are +-(gamma + delta * i), where
    // gamma =          sqrt(( a + sqrt(a * a + b * b)) / 2)
    // delta = sgn(b) * sqrt((-a + sqrt(a * a + b * b)) / 2)
    double r = sqrt(a * a + b * b);
    if (b > 0)
    {
        return ComplexNumber(sqrt((a + r) / 2), sqrt((-a + r) / 2));
    }
    else if (b < 0)
    {
        return ComplexNumber(sqrt((a + r) / 2), -sqrt((-a + r) / 2));
    }
    else // b = 0
    {
        if (a >= 0)
            return ComplexNumber(sqrt(a), 0);
        else // a < 0
            return ComplexNumber(0, sqrt(-a));
    }
}

ComplexNumber ComplexNumber::Euler(double A, double phi)
{
    // Calculate A * e ^ (i * phi)
    //         = A * cos(phi) * A * i * sin(phi)
    return ComplexNumber(A * cos(phi), A * sin(phi));
}

ComplexVector::ComplexVector(ComplexNumber x, ComplexNumber y, ComplexNumber z)
    : x(x), y(y), z(z)
{
}

ComplexVector ComplexVector::operator+(const ComplexVector &v)
{
    return ComplexVector(x + v.x, y + v.y, z + v.z);
}
