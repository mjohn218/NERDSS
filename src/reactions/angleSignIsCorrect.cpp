#include "reactions/association/association.hpp"

bool angleSignIsCorrect(const Vector& vec1, const Vector& vec2)
{
    /*if their normal points in -z, keep theta, otherwise flip sign                                                                                                                                             
    ADDED: Do not flip the sign if it is PI or Zero.                                                                                                                                                          
    double *test=new double[3];
    crossproduct(sigma, normal, test);
    if(test[2]>0  and abs(theta) > 1E-12 and (M_PI - abs(theta)) > 1E-12)//positive z, flip theta.                                                                                                              
	theta=-theta;
    */
    if (std::abs(vec1.z) > 1E-12 && std::abs(vec2.z) > 1E-12) {
        Vector newVec1 { vec1.x, 0, vec1.z };
        Vector newVec2 { vec2.x, 0, vec2.z };
        return newVec1.cross(newVec2).y < 0;
    } else {
        Vector newVec1 { vec1.x, vec1.y, 0 };
        Vector newVec2 { vec2.x, vec2.y, 0 };
        return newVec1.cross(newVec2).z < 0;
    }
}
