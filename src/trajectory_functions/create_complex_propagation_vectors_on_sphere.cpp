#include "math/rand_gsl.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// only works for complex on sphere surface
Coord create_complex_propagation_vectors_on_sphere(const Parameters& params, Complex& targCom)
{
    // TRACE();
    Coord trajTrans;

    double dx = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    double dy = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    double dl = sqrt(dx * dx + dy * dy); // propagation length
    double dangle = acos(dx / dl);
    if (dy < 0.0) {
        dangle = 2.0 * M_PI - dangle;
    } // propagation direction
    double rotangle = dangle - M_PI / 2.0;
    Coord COM = targCom.comCoord;
    Coord COMsphere = find_spherical_coords(COM);
    double dtheta = dl / COMsphere.z;
    Coord COMnewtmp = Coord { COMsphere.x - dtheta, COMsphere.y, COMsphere.z };
    COMnewtmp = find_cardesian_coords(COMnewtmp);
    // define the inner-coords-set
    Vector i = Vector { COM.x, COM.y, COM.z };
    i.normalize();
    Vector temp = Vector { 0.0, 0.0, COM.get_magnitude() };
    temp.normalize();
    Vector j = temp.cross(i);
    j.normalize();
    Vector k = i.cross(j);
    k.normalize();
    //std::vector<Vector> crdset;
    std::array<double, 9> crdset;
    crdset[0] = i.x;
    crdset[1] = i.y;
    crdset[2] = i.z;
    crdset[3] = j.x;
    crdset[4] = j.y;
    crdset[5] = j.z;
    crdset[6] = k.x;
    crdset[7] = k.y;
    crdset[8] = k.z;
    Coord COMnew = rotate_on_sphere(COMnewtmp, COM, crdset, rotangle);

    trajTrans.x = COMnew.x - targCom.comCoord.x;
    trajTrans.y = COMnew.y - targCom.comCoord.y;
    trajTrans.z = COMnew.z - targCom.comCoord.z;
    return trajTrans;
}