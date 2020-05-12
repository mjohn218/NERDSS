/*! \file class_Vector.cpp
* \ingroup Associate
 * Created on 5/20/18 by Matthew Varga
 * Purpose:
 * Notes:
 */

#include "classes/class_Vector.hpp"

#include <cmath>

/* CONSTRUCTORS */
Vector::Vector(double x, double y)
{
    this->x = x;
    this->y = y;
}

Vector::Vector(const double& x, const double& y, const double& z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::Vector(const Coord& coordEnd, const Coord& coordStart)
{
    /// constructor which takes in two points and outputs a vector with the start point at the origin
    this->x = coordEnd.x - coordStart.x;
    this->y = coordEnd.y - coordStart.y;
    this->z = coordEnd.z - coordStart.z;
}

Vector::Vector(Coord coord)
{
    /// constructor where the start point is implied to be the origin
    x = coord.x;
    y = coord.y;
    z = coord.z;
}

Vector::Vector(std::array<double, 3>& arr)
    : Coord(arr)
{
    x = arr[0];
    y = arr[1];
    z = arr[2];
}

Vector::Vector(std::vector<double> arr)
{
    if (arr.size() != 3) {
        std::cerr << "ERROR: Attempting to create Vector from an array with size > 3.";
        exit(1);
    }

    x = arr[0];
    y = arr[1];
    z = arr[2];

    //    magnitude = sqrt(x*x + y*y + z*z);
}

/* OPERATORS */
std::ostream& operator<<(std::ostream& os, Vector& vec)
{
    return os << '[' << vec.x << "i + " << vec.y << "j + " << vec.z << "k]";
}
Coord operator+(const Vector& vec, const Coord& crd) { return { vec.x + crd.x, vec.y + crd.y, vec.z + crd.z }; }

/* MEMBER FUNCTIONS */
Vector Vector::cross(const Vector& vec) const
{
    double u0 = this->y * vec.z - this->z * vec.y;
    double u1 = this->z * vec.x - this->x * vec.z;
    double u2 = this->x * vec.y - this->y * vec.x;

    Vector test { u0, u1, u2 };
    test.normalize();
    return test;
}

double Vector::dot(const Vector& vec) const { return ((this->x * vec.x) + (this->y * vec.y) + (this->z * vec.z)); }

double Vector::dot_theta(const Vector& vec) const
{
    double dp { this->dot(vec) };
    double cTheta = dp / (this->magnitude * vec.magnitude);

    if (std::abs(cTheta - 1) < 1E-12)
        cTheta = 1;
    if (std::abs(-1 - cTheta) < 1E-12)
        cTheta = -1;

    if (this->magnitude < 1E-8 || vec.magnitude < 1E-8) {
        std::cout << "WARNING: Attempted dot product with vector of magnitude 0.\n";
        return 0.0;
    } else
        return acos(cTheta);
}

Vector Vector::vector_projection(Vector normal)
{
    double coefficient { this->dot(normal) / normal.dot(normal) };
    Vector sTerm { normal * coefficient };
    sTerm.calc_magnitude();
    return { *this - sTerm };
}

void Vector::normalize()
{
    calc_magnitude();
    if (magnitude == 0) {
        //dividing by zero results in NANs
        magnitude = 1;
        //    std::cout <<"WARNING: NORMALIZING A VECTOR OF ZERO LENGTH "<<std::endl;
    }
    *this = *this /= magnitude;
    calc_magnitude();
}

void Vector::calc_magnitude()
{
    this->magnitude = sqrt(x * x + y * y + z * z);
}
