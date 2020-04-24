/*! \file class_vector.hpp
 * \ingroup Associate
 * Created on 5/20/18 by Matthew Varga
 * Purpose: Vector class for association
 * Notes:
 */
#pragma once

#include "class_Coord.hpp"

/*! \ingroup Associate
 * \brief Holds a vector with the origin as the start point
 *
 * A 2D vector will always be an x and y coordinate
 */
struct Vector : public Coord {
    double magnitude { 0.0 };

    // operators
    Vector operator*(const double val) { return { x * val, y * val, z * val }; }
    Vector operator-(const double val) { return { x - val, y - val, z - val }; };
    Vector operator-(const Vector& vec) { return { this->x - vec.x, this->y - vec.y, this->z - vec.z }; }
    Vector operator-() { return { -x, -y, -z }; }
    Vector operator/=(double val) { return { x / val, y / val, z / val }; }
    Vector operator/(double val) { return { x / val, y / val, z / val }; }
    Vector operator+(const Coord& crd) { return { this->x + crd.x, this->y + crd.y, this->z + crd.z }; }
    friend std::ostream& operator<<(std::ostream& os, Vector& vec);
    friend Coord operator+(const Vector& vec, const Coord& crd);

    //void operator=(const Vector v)
    //{
    //    this->x = v.x;
    //    this->y = v.y;
    //    this->z = v.z;
    //}

    /*!
     * \brief Calculates the magnitude of the Vector. Note that this is no longer done when the Vector is constructed.
     */
    void calc_magnitude();

    /*!
     * \brief Normalizes the vector to unity.
     */
    void normalize();

    /*!
     * \brief Returns a Vectorwhich is the crossproduct of the Vector with another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    Vector cross(const Vector& vec) const;

    /*!
     * \brief Returns the dot product of the Vector with another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    double dot(const Vector& vec) const; //!< returns the dot product of the vector and another

    /*!
     * \brief Returns the angle between the Vector and another Vector.
     * \params[in] this Vector 1.
     * \params[in] vec Vector 2.
     */
    double dot_theta(const Vector& vec) const; //!< returns the angle between two vectors

    /*!
     * \brief returns the projection of the vector onto the normal.
     *
     * \$f v2 = v1 - \frac{v1 \dot v2}{v2 \dot v2} \times v2 \$f
     */
    Vector vector_projection(Vector normal);

    Vector() = default;
    Vector(double x, double y);
    Vector(const double& x, const double& y, const double& z);
    explicit Vector(std::array<double, 3>& arr);
    explicit Vector(std::vector<double> arr);
    explicit Vector(Coord coord);
    Vector(const Coord& coordEnd, const Coord& coordStart);
};
