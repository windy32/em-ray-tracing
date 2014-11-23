#include "Matrix.h"

Matrix::Matrix(
    double m11, double m12, double m13, 
    double m21, double m22, double m23,
    double m31, double m32, double m33)
{
    this->m11 = m11;
    this->m12 = m12;
    this->m13 = m13;

    this->m21 = m21;
    this->m22 = m22;
    this->m23 = m23;

    this->m31 = m31;
    this->m32 = m32;
    this->m33 = m33;
}

double Matrix::det() const
{
    return 
        m11 * m22 * m33 + 
        m12 * m23 * m31 + 
        m13 * m21 * m32 - 
        m13 * m22 * m31 - 
        m11 * m23 * m32 - 
        m12 * m21 * m33;
}

Matrix Matrix::inverse() const
{
    double d = det();

    double t11 = m22 * m33 - m23 * m32;
    double t12 = m13 * m32 - m12 * m33;
    double t13 = m12 * m23 - m13 * m22;

    double t21 = m23 * m31 - m21 * m33;
    double t22 = m11 * m33 - m13 * m31;
    double t23 = m13 * m21 - m11 * m23;

    double t31 = m21 * m32 - m22 * m31;
    double t32 = m12 * m31 - m11 * m32;
    double t33 = m11 * m22 - m12 * m21;

    return Matrix(
        t11 / d, t12 / d, t13 / d,
        t21 / d, t22 / d, t23 / d,
        t31 / d, t32 / d, t33 / d);
}

Point Matrix::operator*(const Point &p) const
{
    Point result;
    result.x = m11 * p.x + m12 * p.y + m13 * p.z;
    result.y = m21 * p.x + m22 * p.y + m23 * p.z;
    result.z = m31 * p.x + m32 * p.y + m33 * p.z;
    return result;
}

Vector Matrix::operator*(const Vector &v) const
{
    Vector result;
    result.x = m11 * v.x + m12 * v.y + m13 * v.z;
    result.y = m21 * v.x + m22 * v.y + m23 * v.z;
    result.z = m31 * v.x + m32 * v.y + m33 * v.z;
    return result;
}

ComplexVector Matrix::operator*(const ComplexVector &v) const
{
    return ComplexVector(
        ComplexNumber(v.x * m11 + v.y * m12 + v.z * m13),
        ComplexNumber(v.x * m21 + v.y * m22 + v.z * m23),
        ComplexNumber(v.x * m31 + v.y * m32 + v.z * m33));
}
