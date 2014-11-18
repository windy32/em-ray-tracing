#include "Vector.h"
#include <math.h>

Vector::Vector(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::Vector(const Point &start, const Point &end)
{
    x = end.x - start.x;
    y = end.y - start.y;
    z = end.z - start.z;
}

double Vector::length() const
{ 
    return sqrt(x * x + y * y + z * z); 
}

double Vector::sqrLength() const
{ 
    return x * x + y * y + z * z; 
}

double & Vector::operator[](int index)
{
    if (index == 0)
        return x;
    else if (index == 1)
        return y;
    else
        return z;
}

const double & Vector::operator[](int index) const
{
    if (index == 0)
        return x;
    else if (index == 1)
        return y;
    else
        return z;
}

Vector Vector::operator+(const Vector &b) const 
{
    return Vector(x + b.x, y + b.y, z + b.z); 
}

Vector Vector::operator-(const Vector &b) const 
{ 
    return Vector(x - b.x, y - b.y, z - b.z); 
}

Vector Vector::operator*(double b) const 
{ 
    return Vector(x * b, y * b, z * b); 
}

Vector Vector::mult(const Vector &b) const
{ 
    return Vector(x * b.x, y * b.y, z * b.z); 
}

Vector& Vector::norm() 
{ 
    return *this = *this * (1 / sqrt(x * x + y * y + z * z)); 
}

double Vector::dot(const Vector &b) const 
{
    return x * b.x + y * b.y + z * b.z; 
}

double Vector::dot(const Point &p) const
{
    return x * p.x + y * p.y + z * p.z;
}

Vector Vector::cross(const Vector &b)
{
    return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
}

double Vector::angleTo(const Vector &b) const
{
    return acos(dot(b) * (1.0f / (length() * b.length())));
}