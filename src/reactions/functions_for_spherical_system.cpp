/*
 * ### Created on 2/05/2020 by Yiben Fu
 * ### Purpose
 * ***
 * all functions that are used in spherical systems
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"

#include <array>
#include <math.h>

double radius(Coord mol)
{
    return sqrt(mol.x * mol.x + mol.y * mol.y + mol.z * mol.z);
}

Coord find_spherical_coords(Coord mol) // mol: cardesian coords, output spherical coords
{
    Coord angles;
    double R = radius(mol);
    if (mol.z == R) {
        angles.x = 0.0;
        angles.y = 0.0;
        angles.z = R;
    } else if (mol.z == -R) {
        angles.x = -M_PI;
        angles.y = 0.0;
        angles.z = R;
    } else {
        double theta = acos(mol.z / R); // acos, gives angle range [0,pi]
        double phi = acos(mol.x / (R * sin(theta)));
        if (std::isnan(phi)) {
            phi = 0.0;
        }
        if (mol.y < 0)
            phi = 2.0 * M_PI - phi;
        angles.x = theta;
        angles.y = phi;
        angles.z = R;
    }
    return angles;
}

Coord find_cardesian_coords(Coord mol) // mol: spherical coords, output cardesian coords
{
    Coord xyz;
    double theta = mol.x;
    double phi = mol.y;
    double R = mol.z;
    xyz.x = R * sin(theta) * cos(phi);
    xyz.y = R * sin(theta) * sin(phi);
    xyz.z = R * cos(theta);
    return xyz;
}

double theta_plus(double theta1, double theta2) // sum of two theta
{
    double sum = theta1 + theta2;
    if (sum > M_PI) {
        sum = 2.0 * M_PI - sum;
    } else if (sum < 0.0) {
        sum = -sum;
    }
    return sum;
}

double phi_plus(double phi1, double phi2) // sum of two phi
{
    double sum = phi1 + phi2;
    if (sum > 2.0 * M_PI) {
        sum = sum - 2.0 * M_PI;
    } else if (sum < 0.0) {
        sum = 2.0 * M_PI + sum;
    }
    return sum;
}

Coord angle_plus(Coord angle1, Coord angle2)
{
    Coord sum;
    double phi = phi_plus(angle1.y, angle2.y);
    double theta = angle1.x + angle2.x;
    if (theta > M_PI) {
        theta = 2.0 * M_PI - theta;
        phi = phi_plus(phi, M_PI);
    } else if (theta < 0.0) {
        theta = -theta;
        phi = phi_plus(phi, -M_PI);
    } else if (theta == 0.0 || theta == M_PI) {
        phi = 0.0;
    }
    sum.x = theta;
    sum.y = phi;
    sum.z = angle1.z;
    return sum;
}

Coord find_position_after_association(double arc1, Coord Iface1, Coord Iface2, double arc_total, double bindRadius)
{
    double x1 = Iface1.x;
    double y1 = Iface1.y;
    double z1 = Iface1.z;
    double x2 = Iface2.x;
    double y2 = Iface2.y;
    double z2 = Iface2.z;

    double R = radius(Iface1);
    // define the unit vector of the plane, Origin-Iface1-Iface2;
    double nx = 1.0;
    double ny = (x2 * z1 - x1 * z2) / (y1 * z2 - y2 * z1);
    double nz = (x2 * y1 - x1 * y2) / (z1 * y2 - z2 * y1);
    double mag = sqrt(nx * nx + ny * ny + nz * nz);
    nx = nx / mag;
    ny = ny / mag;
    nz = nz / mag;
    /////////////////////////////////////////////////////////
    /* calculate the new position of Iface1.*/
    arc1 = std::abs(arc1);
    double a1 = nz * x1 - nx * z1;
    double a11 = R * R * nz * cos(arc1 / R);
    double a2 = ny * x1 - nx * y1;
    double a22 = R * R * ny * cos(arc1 / R);
    double a3 = ny * z1 - nz * y1;
    // function: A*x^2 + B*x + C = 0;
    double A = a1 * a1 + a2 * a2 + a3 * a3;
    double B = -2.0 * (a1 * a11 + a2 * a22);
    double C = a11 * a11 + a22 * a22 - a3 * a3 * R * R;
    double delta = B * B - 4.0 * A * C;
    if (delta < 0.0)
        delta = 0.0;
    // solution1
    double X1 = 0.5 / A * (-B + sqrt(delta));
    double Y1 = (a1 * X1 - a11) / a3;
    double Z1 = -(a2 * X1 - a22) / a3;
    double distance1 = sqrt(pow(X1 - x2, 2.0) + pow(Y1 - y2, 2.0) + pow(Z1 - z2, 2.0));
    // solution2
    double X2 = 0.5 / A * (-B - sqrt(delta));
    double Y2 = (a1 * X2 - a11) / a3;
    double Z2 = -(a2 * X2 - a22) / a3;
    double distance2 = sqrt(pow(X2 - x2, 2.0) + pow(Y2 - y2, 2.0) + pow(Z2 - z2, 2.0));
    // determine which solution is correct
    Coord new_position1;
    if (bindRadius < arc_total) {
        if (distance1 < distance2) {
            new_position1.x = X1;
            new_position1.y = Y1;
            new_position1.z = Z1;
        } else {
            new_position1.x = X2;
            new_position1.y = Y2;
            new_position1.z = Z2;
        }
    } else {
        if (distance1 > distance2) {
            new_position1.x = X1;
            new_position1.y = Y1;
            new_position1.z = Z1;
        } else {
            new_position1.x = X2;
            new_position1.y = Y2;
            new_position1.z = Z2;
        }
    }
    if (std::isnan(new_position1.x) || std::isnan(new_position1.y)) {
        std::cout << "WRONG: non position is generated in 'find_position_after_association'...EXIT! " << std::endl;
        exit(1);
    }
    return new_position1;
}

/*dtheta is the polar angle change, not the solid angle
  COM and targ has already been moved along the azimuth by dphi.
*/
// COM COMnew are cardeseian coords
std::array<double, 9> inner_coord_set(Coord com, Coord comnew)
{
    std::array<double, 9> crdset;
    crdset[0] = 0;
    crdset[1] = 0;
    crdset[2] = 0;
    crdset[3] = 0;
    crdset[4] = 0;
    crdset[5] = 0;
    crdset[6] = 0;
    crdset[7] = 0;
    crdset[8] = 0;
    Vector i, j, k, v;
    if (sqrt(pow(com.x - comnew.x, 2.0) + pow(com.y - comnew.y, 2.0) + pow(com.z - comnew.z, 2.0)) < 1E-8) {
        //std::cout<<"During the tranlsation on sphere, the complex doesn't move, then use the default inner_coord_set"<<std::endl;
        i = Vector { com.x, com.y, com.z };
        i.normalize();
        v = Vector { 0.0, 0.0, 1.0 };
        if (std::abs(std::abs(com.z) - com.get_magnitude()) < 1E-8) {
            v = Vector { -1.0, 0.0, 0.0 };
        }
        j = v.cross(i);
        k = i.cross(j);
    } else {
        i = Vector { com.x, com.y, com.z };
        i.normalize();
        v = Vector { comnew.x, comnew.y, comnew.z };
        v.normalize();
        k = i.cross(v);
        j = k.cross(i);
    }
    i.normalize();
    j.normalize();
    k.normalize();

    crdset[0] = i.x;
    crdset[1] = i.y;
    crdset[2] = i.z;
    crdset[3] = j.x;
    crdset[4] = j.y;
    crdset[5] = j.z;
    crdset[6] = k.x;
    crdset[7] = k.y;
    crdset[8] = k.z;

    return crdset;
}
std::array<double, 9> inner_coord_set_new(Coord com, Coord comnew)
{
    std::array<double, 9> crdsetnew;
    crdsetnew[0] = 0;
    crdsetnew[1] = 0;
    crdsetnew[2] = 0;
    crdsetnew[3] = 0;
    crdsetnew[4] = 0;
    crdsetnew[5] = 0;
    crdsetnew[6] = 0;
    crdsetnew[7] = 0;
    crdsetnew[8] = 0;
    Vector i, j, k, v;
    if (sqrt(pow(com.x - comnew.x, 2.0) + pow(com.y - comnew.y, 2.0) + pow(com.z - comnew.z, 2.0)) < 1E-8) {
        //std::cout<<"During the tranlsation on sphere, the complex doesn't move, then use the default inner_coord_set_new"<<std::endl;
        i = Vector { comnew.x, comnew.y, comnew.z };
        i.normalize();
        v = Vector { 0.0, 0.0, 1.0 };
        if (std::abs(std::abs(comnew.z) - comnew.get_magnitude()) < 1E-8) {
            v = Vector { -1.0, 0.0, 0.0 };
        }
        j = v.cross(i);
        k = i.cross(j);
    } else {
        Coord dist = comnew - com;
        double l = dist.get_magnitude();
        double R = com.get_magnitude();
        double lnew = l + l * R * R / (R * R - l * l);
        Coord dREFnew = (lnew / l) * dist;
        Coord REFnew = com + dREFnew;
        REFnew = (R / REFnew.get_magnitude()) * REFnew;
        com = comnew;
        comnew = REFnew;

        i = Vector { com.x, com.y, com.z };
        v = Vector { comnew.x, comnew.y, comnew.z };
        k = i.cross(v);
        j = k.cross(i);
    }
    i.normalize();
    j.normalize();
    k.normalize();

    crdsetnew[0] = i.x;
    crdsetnew[1] = i.y;
    crdsetnew[2] = i.z;
    crdsetnew[3] = j.x;
    crdsetnew[4] = j.y;
    crdsetnew[5] = j.z;
    crdsetnew[6] = k.x;
    crdsetnew[7] = k.y;
    crdsetnew[8] = k.z;

    return crdsetnew;
}
// crdset is the previous one, not the new or updated one
std::array<double, 3> calculate_inner_coord_coefficients(Coord TARG, Coord COM, std::array<double, 9> crdset)
{
    std::array<double, 3> M {};
    M[0] = 0.0;
    M[1] = 0.0;
    M[2] = 0.0;

    Vector targ = Vector(TARG - COM);
    targ.calc_magnitude();
    if (targ.magnitude < 1E-8) { // targ is as com
        return M;
    }

    // get the inner_coords_set
    Vector i = Vector { crdset[0], crdset[1], crdset[2] };
    Vector j = Vector { crdset[3], crdset[4], crdset[5] };
    Vector k = Vector { crdset[6], crdset[7], crdset[8] };

    // targ = alpha*i + beta*j + gama*k;
    double alpha, beta, gama;
    // check whether targ vector is perpenticular to any of the coords set_memProtein_sphere
    Vector Targ = targ;
    Targ.normalize();
    if (std::abs(Targ.dot(i)) < 1E-8) { // targ is perpenticular to i
        alpha = 0.0;
        if (std::abs(Targ.dot(j)) < 1E-8) { // targ is also perpenticular to j
            beta = 0.0;
            gama = targ.magnitude;
            if (Targ.dot(k) < 0.0) {
                gama = -gama;
            }
        } else if (std::abs(Targ.dot(k)) < 1E-8) { // targ is also perpenticular to k
            gama = 0.0;
            beta = targ.magnitude;
            if (Targ.dot(j) < 0.0) {
                beta = -beta;
            }
        } else {
            beta = (targ.x * k.y - targ.y * k.x) / (j.x * k.y - j.y * k.x);
            gama = (targ.x * j.y - targ.y * j.x) / (k.x * j.y - k.y * j.x);
        }
    } else if (std::abs(Targ.dot(j)) < 1E-8) { // targ is verticle to j
        beta = 0.0;
        if (std::abs(Targ.dot(i)) < 1E-8) { // verticle to i
            alpha = 0.0;
            gama = targ.magnitude;
            if (Targ.dot(k) < 0.0) {
                gama = -gama;
            }
        } else if (std::abs(Targ.dot(k)) < 1E-8) { // verticle to k
            gama = 0.0;
            alpha = targ.magnitude;
            if (Targ.dot(i) < 0.0) {
                alpha = -alpha;
            }
        } else {
            alpha = (targ.x * k.y - targ.y * k.x) / (i.x * k.y - i.y * k.x);
            gama = (targ.x * i.y - targ.y * i.x) / (k.x * i.y - k.y * i.x);
        }
    } else if (std::abs(Targ.dot(k)) < 1E-8) { // targ is verticle to k
        gama = 0.0;
        if (std::abs(Targ.dot(i)) < 1E-8) { // verticle to i
            alpha = 0.0;
            beta = targ.magnitude;
            if (Targ.dot(j) < 0.0) {
                beta = -beta;
            }
        } else if (std::abs(Targ.dot(j)) < 1E-8) { // verticle to j
            beta = 0.0;
            alpha = targ.magnitude;
            if (Targ.dot(i) < 0.0) {
                alpha = -alpha;
            }
        } else {
            alpha = (targ.x * j.y - targ.y * j.x) / (i.x * j.y - i.y * j.x);
            beta = (targ.x * i.y - targ.y * i.x) / (j.x * i.y - j.y * i.x);
        }
    } else {
        double n1 = k.y * targ.x - k.x * targ.y;
        double n2 = i.x * k.y - i.y * k.x;
        double n3 = j.x * k.y - j.y * k.x;
        double n4 = k.z * targ.x - k.x * targ.z;
        double n5 = i.x * k.z - i.z * k.x;
        double n6 = j.x * k.z - j.z * k.x;
        alpha = (n1 * n6 - n4 * n3) / (n2 * n6 - n5 * n3);
        beta = (n1 * n6 - n2 * n6 * alpha) / (n3 * n6);
        gama = (targ.x - alpha * i.x - beta * j.x) / k.x;
    }
    M[0] = alpha;
    M[1] = beta;
    M[2] = gama;
    return M;
}

// input and output are cardesian coords
Coord translate_on_sphere(Coord targ, Coord COM, Coord COMnew, std::array<double, 9> crdset, std::array<double, 9> crdsetnew)
{
    Coord dcom = COMnew - COM;
    if (dcom.get_magnitude() < 1E-8) { // no translation on sphere
        return targ;
    }
    // calculate the inner-coords-set efficient
    std::array<double, 3> M = calculate_inner_coord_coefficients(targ, COM, crdset);
    double alpha = M[0];
    double beta = M[1];
    double gama = M[2];
    Vector i = Vector { crdsetnew[0], crdsetnew[1], crdsetnew[2] };
    Vector j = Vector { crdsetnew[3], crdsetnew[4], crdsetnew[5] };
    Vector k = Vector { crdsetnew[6], crdsetnew[7], crdsetnew[8] };
    Coord targnew = Coord { alpha * i + beta * j + gama * k };
    targnew = targnew + COMnew;
    return targnew;
}

// input and output are cardesian coords
Coord rotate_on_sphere(Coord Targ, Coord COM, std::array<double, 9> crdset, double dangle)
{
    Coord targnew;

    // targ will rotate along O-COM line, i.e. i axis. so its projection along i is not change
    Vector i = Vector { crdset[0], crdset[1], crdset[2] };
    Vector j = Vector { crdset[3], crdset[4], crdset[5] };
    Vector k = Vector { crdset[6], crdset[7], crdset[8] };
    Vector targ = Vector { Targ - COM };

    Vector targi = i * targ.dot(i);
    Vector targjk = targ - targi;
    targi.calc_magnitude();
    targjk.calc_magnitude();
    // if targ is on the line of O-COM, then no need to rotate
    if (targjk.magnitude < 1E-8 || std::abs(targi.magnitude - 1.0) < 1E-8) {
        targnew = Targ;
    } else {
        double phi = acos(targjk.dot(j) / targjk.magnitude);
        if (targjk.dot(k) < 0.0) {
            phi = 2.0 * M_PI - phi;
        }
        phi = phi + dangle;
        targjk = j * (targjk.magnitude * cos(phi)) + k * (targjk.magnitude * sin(phi));
        targ = targi + targjk;
        targnew = Coord { targ.x + COM.x, targ.y + COM.y, targ.z + COM.z };
    }
    if (std::isnan(targnew.x)) {
        if (targjk.magnitude < 1E-8 || std::abs(targi.magnitude - 1.0) < 1E-8) {
            targnew = Coord { Targ.x, Targ.y, Targ.z };
        } else {
            std::cout << "WRONG! NON is generated after the rotation on sphere! EXIT..." << std::endl;
            exit(1);
        }
    }
    return targnew;
}

double calc_bindRadius2D(double bindRadius, Coord iFace)
{
    double R;
    R = radius(iFace);
    double bindRadius2D;
    bindRadius2D = R * 2 * asin((0.5 * bindRadius) / R);
    return bindRadius2D;
}

void set_memProtein_sphere(Complex reactCom, Molecule& memProtein, std::vector<Molecule> moleculeList, const Membrane membraneObject)
{
    //if (membraneObject.implicitLipid == false){ //for explicit lipid model; lipid is a member of reactCom
    double r = 0.0;
    for (auto mol : reactCom.memberList) {
        if ((moleculeList[mol].isLipid == true || moleculeList[mol].isImplicitLipid == true) && moleculeList[mol].tmpComCoord.get_magnitude() > r) {
            memProtein = moleculeList[mol];
            r = moleculeList[mol].tmpComCoord.get_magnitude();
        }
    }
    memProtein.comCoord = memProtein.tmpComCoord;
    for (int i = 0; i < memProtein.interfaceList.size(); i++) { // here memProtein is an Lipid, has only one interface
        Coord ifaceToCom = memProtein.tmpICoords[i] - memProtein.tmpComCoord;
        double bond = ifaceToCom.get_magnitude();
        double R = memProtein.comCoord.get_magnitude();
        memProtein.interfaceList[i].coord = (R - bond) / R * memProtein.comCoord;
    }

    // for implicit lipid model, on 2D->2D case, reactCom has no implicitlipid member.
    if (memProtein.isLipid == false && memProtein.isImplicitLipid == false && membraneObject.implicitLipid == true) {
        Coord iface;
        Coord com;
        r = 0.0;
        for (auto mol : reactCom.memberList) {
            for (int i = 0; i < moleculeList[mol].interfaceList.size(); i++) {
                if (moleculeList[mol].interfaceList[i].isBound == true) {
                    int index = moleculeList[mol].interfaceList[i].interaction.partnerIndex;
                    if (moleculeList[index].isImplicitLipid == true && moleculeList[mol].tmpICoords[i].get_magnitude() > r) {
                        memProtein = moleculeList[index];
                        com = moleculeList[mol].tmpComCoord;
                        iface = moleculeList[mol].tmpICoords[i];
                        r = moleculeList[mol].tmpICoords[i].get_magnitude();
                    }
                }
            }
        }
        memProtein.comCoord = iface; // here targ is an ImplicitLipid, has only one interface
        Coord ifaceToCom = iface - com;
        double bond = ifaceToCom.get_magnitude();
        double R = memProtein.comCoord.get_magnitude();
        memProtein.interfaceList[0].coord = (R - bond) / R * memProtein.comCoord;
    }
    if (memProtein.isLipid == false && memProtein.isImplicitLipid == false) {
        std::cout << "WRONG: failed to create memProtein, in the step to adjust complex's orientation on sphere. Exit..." << std::endl;
        exit(1);
    }
    //memProtein.comCoord =  (memProtein.comCoord.get_magnitude() + 0.1) / memProtein.comCoord.get_magnitude()  * memProtein.comCoord;
    //memProtein.interfaceList[0].coord = (memProtein.interfaceList[0].coord.get_magnitude() + 0.1)/memProtein.interfaceList[0].coord.get_magnitude() *  memProtein.interfaceList[0].coord;
}
void find_Lipid_sphere(Complex reactCom, Molecule& Lipid, std::vector<Molecule> moleculeList, const Membrane membraneObject)
{
    //if (membraneObject.implicitLipid == false){ //for explicit lipid model; lipid is a member of reactCom
    double r = 0.0;
    for (auto mol : reactCom.memberList) {
        if ((moleculeList[mol].isLipid == true || moleculeList[mol].isImplicitLipid == true) && moleculeList[mol].tmpComCoord.get_magnitude() > r) {
            Lipid = moleculeList[mol];
            r = moleculeList[mol].tmpComCoord.get_magnitude();
        }
    }
    // for implicit lipid model, on 2D->2D case, reactCom has no implicitlipid member.
    if (Lipid.isLipid == false && Lipid.isImplicitLipid == false && membraneObject.implicitLipid == true) {
        r = 0.0;
        Coord iface;
        Coord com;
        for (auto mol : reactCom.memberList) {
            for (int i = 0; i < moleculeList[mol].interfaceList.size(); i++) {
                if (moleculeList[mol].interfaceList[i].isBound == true) {
                    int index = moleculeList[mol].interfaceList[i].interaction.partnerIndex;
                    if (moleculeList[index].isImplicitLipid == true && moleculeList[mol].tmpICoords[i].get_magnitude() > r) {
                        Lipid = moleculeList[index];
                        com = moleculeList[mol].tmpComCoord;
                        iface = moleculeList[mol].tmpICoords[i];
                        r = moleculeList[mol].tmpICoords[i].get_magnitude();
                    }
                }
            }
        }
        Lipid.set_tmp_association_coords();
        Lipid.comCoord = iface; // here targ is an ImplicitLipid, has only one interface
        Lipid.tmpComCoord = iface;
        Lipid.interfaceList[0].coord = com;
        Lipid.tmpICoords[0] = com;
    }

    if (Lipid.isLipid == false && Lipid.isImplicitLipid == false) {
        std::cout << "WRONG: failed to create memProtein, in the step to adjust complex's orientation on sphere. Exit..." << std::endl;
        exit(1);
    }

    //Lipid.comCoord = ( Lipid.comCoord.get_magnitude()+ 0.1)/Lipid.comCoord.get_magnitude() * Lipid.comCoord;
    //Lipid.tmpComCoord = ( Lipid.tmpComCoord.get_magnitude()+ 0.1)/Lipid.tmpComCoord.get_magnitude() * Lipid.tmpComCoord;
    //Lipid.interfaceList[0].coord = ( Lipid.interfaceList[0].coord.get_magnitude()+ 0.1)/Lipid.interfaceList[0].coord.get_magnitude() * Lipid.interfaceList[0].coord;
    //Lipid.tmpICoords[0] = ( Lipid.tmpICoords[0].get_magnitude()+ 0.1)/Lipid.tmpICoords[0].get_magnitude() * Lipid.tmpICoords[0];
}