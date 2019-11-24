/*! \file matrix_functions.cpp

 * ### Created on 11/2/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "math/matrix.hpp"
#include <cmath>

Vector matrix_rotate(Vector& vec, std::array<double, 9>& M)
{
    return { M[0] * vec.x + M[1] * vec.y + M[2] * vec.z,
             M[3] * vec.x + M[4] * vec.y + M[5] * vec.z,
             M[6] * vec.x + M[7] * vec.y + M[8] * vec.z };
}

std::array<double, 9> create_euler_rotation_matrix(double x, double y, double z)
{
    double sx{ sin(x) };
    double cx{ cos(x) };
    double sy{ sin(y) };
    double cy{ cos(y) };
    double sz{ sin(z) };
    double cz{ cos(z) };

    std::array<double, 9> M{};
    M[0] = cz * cy;
    M[1] = cz * sx * sy - sz * cx;
    M[2] = cz * sy * cx + sz * sx;
    M[3] = sz * cy;
    M[4] = sz * sx * sy + cz * cx;
    M[5] = sz * sy * cx - cz * sx;
    M[6] = -sy;
    M[7] = cy * sx;
    M[8] = cy * cx;

    return M;
}

std::array<double, 9> create_euler_rotation_matrix(const Coord& angles)
{
    double sx{ sin(angles.x) };
    double cx{ cos(angles.x) };
    double sy{ sin(angles.y) };
    double cy{ cos(angles.y) };
    double sz{ sin(angles.z) };
    double cz{ cos(angles.z) };

    std::array<double, 9> M{};
    M[0] = cz * cy;
    M[1] = cz * sx * sy - sz * cx;
    M[2] = cz * sy * cx + sz * sx;
    M[3] = sz * cy;
    M[4] = sz * sx * sy + cz * cx;
    M[5] = sz * sy * cx - cz * sx;
    M[6] = -sy;
    M[7] = cy * sx;
    M[8] = cy * cx;

    return M;
}
