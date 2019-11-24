/*! \file class_Quat.hpp
 * \brief Header file for Quat class.
 */


#pragma once
#include "classes/class_Vector.hpp"
#include <iostream>

struct Quat {
    double w{ 0 };
    double x{ 0 };
    double y{ 0 };
    double z{ 0 };

    Quat operator*(const Quat& q);
    friend std::ostream& operator<<(std::ostream& os, const Quat& q);
    double norm();
    double mag();
    Quat scale(double scal);
    Quat inverse();

    /*!
     * \brief Returns a Quat with magnitude unity.
     */
    Quat unit();

    /*!
     * \brief Takes the conjugate of the Quat, i.e. \$f Q^* \$f.
     */
    Quat conjugate();

    /*!
     * \brief Performs a vector rotation with a quaternion. See the \ref association page.
     */
    void rotate(Vector& vec);

    Quat() = default;
    Quat(double _w, double _x, double _y, double _z)
        : w(_w)
        , x(_x)
        , y(_y)
        , z(_z)
    {
    }
};
