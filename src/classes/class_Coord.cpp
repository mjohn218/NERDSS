/*! \file class_Coord.cpp

 * Created on 6/1/18 by Matthew Varga
 * Purpose:
 * Notes:
 */

#include "classes/class_Coord.hpp"
#include <classes/class_Parameters.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

// CONSTRUCTORS //
Coord::Coord(std::array<double, 3>& arr)
    : x(arr[0])
    , y(arr[1])
    , z(arr[2])
{
}

Coord::Coord(double x, double y, double z)
    : x(x)
    , y(y)
    , z(z)
{
}

Coord::Coord(std::vector<double> vals)
    : x(vals[0])
    , y(vals[1])
    , z(vals[2])
{
    if (vals.size() > 3) {
        std::cout << "Coordinate cannot have more than 3 points. Exiting." << '\n';
        exit(1);
    }
}

Coord round(const Coord& c) { return { roundv(c.x), roundv(c.y), roundv(c.z) }; }

double roundv(double var)
{
    // if-else is because neg and pos values will round differently
    double val = (int)(var > 0 ? var * 10000 + 0.5 : var * 10000 - 0.5);
    return val / 10000;
}

bool is_co_linear(const Coord& c1, const Coord& c2, const Coord& c3)
{
    bool isCoLinear { false };
    //Heron's formula: S^2=p(p-a)(p-b)(p-c); if the three points co-linear, the S^2 == 0
    double a { 0.0 };
    double b { 0.0 };
    double c { 0.0 };
    a = pow((c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y) + (c1.z - c2.z) * (c1.z - c2.z), 0.5);
    b = pow((c1.x - c3.x) * (c1.x - c3.x) + (c1.y - c3.y) * (c1.y - c3.y) + (c1.z - c3.z) * (c1.z - c3.z), 0.5);
    c = pow((c3.x - c2.x) * (c3.x - c2.x) + (c3.y - c2.y) * (c3.y - c2.y) + (c3.z - c2.z) * (c3.z - c2.z), 0.5);
    double p { (a + b + c) / 2.0 };
    double S { pow(p * (p - a) * (p - b) * (p - c), 0.5) };
    if (std::abs(S - 0.0) < 1E-8) {
        isCoLinear = true;
    }
    return isCoLinear;
}

// OPERATORS //
bool operator==(const Coord& c1, const Coord& c2)
{
    return (bool { roundv(c1.x) == roundv(c2.x) } * bool { roundv(c1.y) == roundv(c2.y) }
        * bool { roundv(c1.z) == roundv(c2.z) });
}

bool operator!=(const Coord& c1, const Coord& c2)
{
    return (bool { roundv(c1.x) != roundv(c2.x) } || bool { roundv(c1.y) != roundv(c2.y) }
        || bool { roundv(c1.z) != roundv(c2.z) });
}

std::ostream& operator<<(std::ostream& os, const Coord& c)
{
    os << std::setprecision(6) << std::setw(12) << c.x << " " << std::setw(12) << c.y << " " << std::setw(12) << c.z;
    return os;
}

Coord operator+(const std::array<double, 3>& arr, const Coord& c)
{
    double x { arr[0] + c.x };
    double y { arr[1] + c.y };
    double z { arr[2] + c.z };
    return { x, y, z };
}

Coord operator+(const Coord& c1, const Coord& c2)
{
    double x { c1.x + c2.x };
    double y { c1.y + c2.y };
    double z { c1.z + c2.z };
    return { x, y, z };
}

Coord operator-(const Coord& c1, const double val) { return Coord { c1.x - val, c1.y - val, c1.z - val }; }

void Coord::operator+=(const Coord& coord)
{
    x += coord.x;
    y += coord.y;
    z += coord.z;
}

void operator/=(Coord& c, double& scal)
{
    c.x = c.x / scal;
    c.y = c.y / scal;
    c.z = c.z / scal;
}

Coord operator+=(Coord& c, const std::array<double, 3>& arr)
{
    double x { c.x + arr[0] };
    double y { c.y + arr[1] };
    double z { c.z + arr[2] };
    return { x, y, z };
}

// MEMBER FUNCTIONS //
void Coord::zero_crds()
{
    /// Just sets all values to zero
    x = 0;
    y = 0;
    z = 0;
}

bool Coord::isOutOfBox(const Membrane& membraneObject)
{
    /*!
     * \brief Checks if the coordinate is outside the waterbox. Only used in reflecting boundary conditions.
     */

    if ((x > (membraneObject.waterBox.x / 2.0)) || (x < -(membraneObject.waterBox.x / 2.0)))
        return true;
    if ((y > (membraneObject.waterBox.y / 2.0)) || (y < -(membraneObject.waterBox.y / 2.0)))
        return true;
    if ((z > (membraneObject.waterBox.z / 2.0)) || (z < -(membraneObject.waterBox.z / 2.0)))
        return true;

    return false;
}

double Coord::get_magnitude()
{
    return sqrt(x * x + y * y + z * z);
}
