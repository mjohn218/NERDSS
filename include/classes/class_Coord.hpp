/*! \file class_coord.hpp

 * Created on 6/1/18 by Matthew Varga
 * Purpose:
 * Notes:
 */
#pragma once

#include "classes/class_Membrane.hpp"
#include "classes/class_Parameters.hpp"

#include <array>
#include <iostream>
#include <vector>

/*! \struct Coord
 * \brief Class to hold xyz coordinates
 */
struct Coord {
public:
    double x { 0 };
    double y { 0 };
    double z { 0 };

    // operator overloads
    void operator+=(const Coord& coord);
    friend bool operator==(const Coord& c1, const Coord& c2);
    friend bool operator!=(const Coord& c1, const Coord& c2);
    friend std::ostream& operator<<(std::ostream& os, const Coord& c);
    friend Coord operator+(const std::array<double, 3>& arr, const Coord& c);
    friend Coord operator+(const Coord& c1, const Coord& c2);
    friend Coord operator-(const Coord& c1, const double val);
    friend Coord operator+=(Coord& c, const std::array<double, 3>& arr);
    friend void operator/=(Coord& c, double& scal);
    void operator=(const std::array<double, 3>& arr)
    {
        this->x = arr[0];
        this->y = arr[1];
        this->z = arr[2];
    }
    //   void operator=(Coord& c)
    //   {
    //       this->x = c.x;
    //       this->y = c.y;
    //       this->z = c.z;
    //   }
    Coord& operator-=(Coord& c)
    {
        this->x -= c.x;
        this->y -= c.y;
        this->z -= c.z;
        return *this;
    }
    Coord operator-(const Coord& coord2) const
    {
        return { this->x - coord2.x, this->y - coord2.y, this->z - coord2.z };
    };

    void zero_crds();
    bool isOutOfBox(const Membrane& membraneObject);

    double get_magnitude();

    Coord() = default;
    Coord(double x, double y, double z);

    // TODO: include this in association 2.0
    explicit Coord(std::array<double, 3>& arr);
    explicit Coord(std::vector<double> vals);
};

Coord round(const Coord& c);
double roundv(double var);

template <typename Scal>
Coord operator*(Scal scal, const Coord& coord)
{
    return { scal * coord.x, scal * coord.y, scal * coord.z };
}

bool is_co_linear(const Coord& c1, const Coord& c2, const Coord& c3);
