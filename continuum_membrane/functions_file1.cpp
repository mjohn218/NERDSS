#include <math.h>
#include <iostream>
#include <armadillo>
#include <omp.h>
//#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#include "functions_file2.cpp"

void element_energy_force_regular(mat dots, Param param, double c0, double& Ebe, mat& F_be, mat& F_s, mat& F_v, rowvec gqcoeff, cube shape_functions) {
    // F_be is the force related to the curvature; F_s is the force related to the area-constraint; F_v is the force related to the volume-constraint.
    // initialize output parameters
    Ebe = 0.0;
    F_be.fill(0);
    F_s.fill(0);
    F_v.fill(0);
    //////////////////////////////////////////////////////////////
    double kc = param.kc;
    double us = param.us;
    double uv = param.uv;
    double S0 = param.S0;
    double V0 = param.V0;
    double S  = param.S;
    double V  = param.V;
    int GaussQuadratureN = param.GaussQuadratureN;
    // Gaussian quadrature, second-order or 3 points.
    mat VWU = setVMU(GaussQuadratureN); 
    //rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    for (int i = 0; i < VWU.n_rows; i++) {
        double ebe = 0.0;
        mat f_be(12,3); f_be.fill(0);
        mat f_cons(12,3); f_cons.fill(0);
        mat f_conv(12,3); f_conv.fill(0);
        //rowvec vwu = VWU.row(i);
        //mat sf(12,7);
        //shapefunctions(vwu,sf);          // 12 shape functions
        mat sf = shape_functions.slice(i);
        // a_1,2,3 covariant vectors; a1,2 contravariant vectors;
        // a_11: a_1 differential to v; a_12: a_1 differential to w;
        rowvec x(3); trans_time(sf.col(0),dots,x);
        rowvec a_1(3); trans_time(sf.col(1),dots,a_1);
        rowvec a_2(3); trans_time(sf.col(2),dots,a_2);
        rowvec a_11(3); trans_time(sf.col(3),dots,a_11);
        rowvec a_22(3); trans_time(sf.col(4),dots,a_22);
        rowvec a_12(3); trans_time(sf.col(5),dots,a_12);
        rowvec a_21(3); trans_time(sf.col(6),dots,a_21);
        rowvec xa = cross(a_1,a_2);
        double sqa = norm(xa,2);
        rowvec xa_1 = cross(a_11,a_2) + cross(a_1,a_21);
        rowvec xa_2 = cross(a_12,a_2) + cross(a_1,a_22);
        mat oneelement1 = xa*strans(xa_1);
        double sqa_1 = 1.0/sqa * oneelement1(0,0);
        mat oneelement2 = xa*strans(xa_2);
        double sqa_2 = 1.0/sqa * oneelement2(0,0);
        rowvec a_3 = xa/sqa;
        rowvec a_31 = 1.0/sqa/sqa *(xa_1*sqa-xa*sqa_1);
        rowvec a_32 = 1.0/sqa/sqa *(xa_2*sqa-xa*sqa_2);
        rowvec d = a_3;
        rowvec d_1 = a_31;
        rowvec d_2 = a_32;
        rowvec a1 = cross(a_2,a_3)/sqa;
        rowvec a2 = cross(a_3,a_1)/sqa;
        rowvec a11 = 1.0/sqa/sqa *( (cross(a_21,a_3)+cross(a_2,a_31))*sqa - cross(a_2,a_3)*sqa_1 );
        rowvec a12 = 1.0/sqa/sqa *( (cross(a_22,a_3)+cross(a_2,a_32))*sqa - cross(a_2,a_3)*sqa_2 );
        rowvec a21 = 1.0/sqa/sqa *( (cross(a_31,a_1)+cross(a_3,a_11))*sqa - cross(a_3,a_1)*sqa_1 );
        rowvec a22 = 1.0/sqa/sqa *( (cross(a_32,a_1)+cross(a_3,a_12))*sqa - cross(a_3,a_1)*sqa_2 );
        mat shu = a1*strans(d_1) + a2*strans(d_2);
        double H_curv = -0.5*shu(0,0);
        ebe = 0.5*kc*sqa*pow(2.0*H_curv-c0,2.0);    // bending energy
        /////////////////////////////////////////////////////////////
        // nodal force
        rowvec n1_be = kc*(2.0*H_curv-c0)*(a1*strans(a1)*d_1+a1*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
        rowvec n2_be = kc*(2.0*H_curv-c0)*(a2*strans(a1)*d_1+a2*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
        rowvec m1_be = -kc*(2.0*H_curv-c0)*a1;
        rowvec m2_be = -kc*(2.0*H_curv-c0)*a2;
        rowvec n1_cons = us*(S-S0)*a1;
        rowvec n2_cons = us*(S-S0)*a2;
        rowvec n1_conv = 1.0/3.0*uv*(V-V0)*(x*strans(d)*a1-x*strans(a1)*d);
        rowvec n2_conv = 1.0/3.0*uv*(V-V0)*(x*strans(d)*a2-x*strans(a2)*d);
        for (int j = 0; j < 12; j++) {
            mat da1 = -sf(j,3)*kron(strans(a1),d) - sf(j,1)*kron(strans(a11),d) - sf(j,1)*kron(strans(a1),d_1) - sf(j,6)*kron(strans(a2),d) - sf(j,2)*kron(strans(a21),d) - sf(j,2)*kron(strans(a2),d_1);
            mat da2 = -sf(j,5)*kron(strans(a1),d) - sf(j,1)*kron(strans(a12),d) - sf(j,1)*kron(strans(a1),d_2) - sf(j,4)*kron(strans(a2),d) - sf(j,2)*kron(strans(a22),d) - sf(j,2)*kron(strans(a2),d_2);
            rowvec tempf_be = n1_be*sf(j,1) + m1_be*da1 + n2_be*sf(j,2) + m2_be*da2;
            f_be.row(j) = tempf_be*sqa;
            rowvec tempf_cons = n1_cons*sf(j,1) + n2_cons*sf(j,2);
            f_cons.row(j) = tempf_cons*sqa;
            rowvec tempf_conv = n1_conv*sf(j,1) + n2_conv*sf(j,2) + 1.0/3.0*uv*(V-V0)*d*sf(j,0);
            f_conv.row(j) = tempf_conv*sqa;
        }
        Ebe = Ebe + 1.0/2.0*gqcoeff(i)*ebe;
        F_be = F_be + 1.0/2.0*gqcoeff(i)*f_be;
        F_s = F_s + 1.0/2.0*gqcoeff(i)*f_cons;
        F_v = F_v + 1.0/2.0*gqcoeff(i)*f_conv;
    }
}

void element_energy_force_irregular(int facenum, Mat<int> A, mat vertex, Param param, double c0, double& E_bending, mat& cunchu1, mat& cunchu2, mat& cunchu3, rowvec gqcoeff, cube shape_functions ) {
    // cunchu1 is the curvature force; cunchu2 is the area-constraint force; cunchu3 is the volume-constraint force.
    // initialize the output parameters
    E_bending = 0.0;
    cunchu1.fill(0.0);
    cunchu2.fill(0.0);
    cunchu3.fill(0.0);
    // five matrix used for subdivision of the irregular patch
    mat M(17,11); M.fill(0.0);
    mat M1(12,17); M1.fill(0.0);
    mat M2(12,17); M2.fill(0.0);
    mat M3(12,17); M3.fill(0.0);
    mat M4(11,17); M4.fill(0.0);
    subdivision_matrix(M, M1, M2, M3, M4);
    int n = 5; // subdivision times
    int i = facenum;
    /////////////////////////////////////////////////////////////////
    // bending energy E_bending and cunchu2_constraint force
    mat ori_dots(11,3);
    for (int k = 0; k < 11; k++) {
        int nodenum = A(i,k+1);
        ori_dots.row(k) = vertex.row(nodenum);
    }
    mat temp = eye(11,11);
    for (int j = 0; j < n; j++) {
        mat newnodes17 = M*ori_dots; // 17 new nodes

        if (j != 0) {
            temp = (M4*M) * temp;
        }
        mat matrix = M*temp;

        mat dots = M1*newnodes17;    // element 1
        mat f1(12,3);  f1.fill(0.0); // f1(12,3)
        mat f2 = f1; 
        mat f3 = f1;
        double ebe1 = 0.0;
        element_energy_force_regular(dots, param, c0, ebe1, f1, f2, f3, gqcoeff, shape_functions);
        mat m1m = strans(M1*matrix);
        cunchu1 = cunchu1 + m1m*f1;
        cunchu2 = cunchu2 + m1m*f2;
        cunchu3 = cunchu3 + m1m*f3;
        E_bending = E_bending + ebe1;

        dots = M2*newnodes17;    // element 2
        f1.fill(0.0);
        f2.fill(0.0);
        f3.fill(0.0);
        double ebe2 = 0.0;
        element_energy_force_regular(dots, param, c0, ebe2, f1, f2, f3, gqcoeff, shape_functions);
        mat m2m = strans(M2*matrix);
        cunchu1 = cunchu1 + m2m*f1;
        cunchu2 = cunchu2 + m2m*f2;
        cunchu3 = cunchu3 + m2m*f3;
        E_bending = E_bending + ebe2;

        dots = M3*newnodes17;    // element 3
        f1.fill(0.0);
        f2.fill(0.0);
        f3.fill(0.0);
        double ebe3 = 0.0;
        element_energy_force_regular(dots, param, c0, ebe3, f1, f2, f3, gqcoeff, shape_functions);
        mat m3m = strans(M3*matrix);
        cunchu1 = cunchu1 + m3m*f1;
        cunchu2 = cunchu2 + m3m*f2;
        cunchu3 = cunchu3 + m3m*f3;
        E_bending = E_bending + ebe3;

        mat dots4 = M4*newnodes17;   // element 4, still irregular patch
        ori_dots = dots4;
    }
}

void energy_force_regularization(mat vertex, mat vertexold, Mat<int> face, Param param, Mat<int> insertionpatch, Row<int> isinsertionpatch, double& E_regularization, mat& Fre, double& E_insertions, rowvec& deformnumbers){
    double k  = param.k;
    double K  = param.K;
    //Fre.fill(0.0);
    E_regularization = 0.0;
    E_insertions = 0.0;
    deformnumbers.fill(0);
    //vec elements_triangle = calculate_all_triangles_area(vertex,face); // this is triangular mesh area, while elementS is limit surface area.
    //double a0 = mean(elements_triangle);
    int deformnumber_shape = 0;
    int deformnumber_area = 0;
    rowvec Ere(face.n_rows); Ere.fill(0.0);
    rowvec Ein(face.n_rows); Ein.fill(0.0);
    //cube fre(vertex.n_rows,3,face.n_rows);
    mat fre(vertex.n_rows,3); fre.fill(0.0);
    #pragma omp parallel for reduction(+:fre)
    for (int i = 0; i < face.n_rows; i++){    
        bool isInsertionPatch = false;
        if ( isinsertionpatch(i) == 1 ){
            isInsertionPatch = true; 
        }
        int node0 = face(i,0); // three nodes of this face element
        int node1 = face(i,1);
        int node2 = face(i,2);
        rowvec vector0 = vertex.row(node0) - vertex.row(node1);  double side0 = norm(vector0,2.0);
        rowvec vector1 = vertex.row(node1) - vertex.row(node2);  double side1 = norm(vector1,2.0);
        rowvec vector2 = vertex.row(node2) - vertex.row(node0);  double side2 = norm(vector2,2.0);
        double s = (side0 + side1 + side2)/2.0;
        double area = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
        double meanside = 1.0/3.0*(side0 + side1+ side2);
        double gama = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
        rowvec vectorold0 = vertexold.row(node0) - vertexold.row(node1);  double sideold0 = norm(vectorold0,2.0);
        rowvec vectorold1 = vertexold.row(node1) - vertexold.row(node2);  double sideold1 = norm(vectorold1,2.0);
        rowvec vectorold2 = vertexold.row(node2) - vertexold.row(node0);  double sideold2 = norm(vectorold2,2.0);
        s = (sideold0 + sideold1 + sideold2)/2.0;
        double areaold = sqrt( s*(s-sideold0)*(s-sideold1)*(s-sideold2) ); //double areaold = S0/face.n_rows;
            
        bool isDeformShape = false;
        if ( gama > param.gama_shape ){
            isDeformShape = true;
        }
        bool isDeformArea = false;
        double a0 = areaold;
        if ( abs(area-a0)/a0 >= param.gama_area ){
            //isDeformArea = true;
        }
            
        if ( isDeformShape == false && isDeformArea == false && isInsertionPatch == false ){
            Ere(i) = k/2.0*(pow(side0-sideold0,2.0) + pow(side1-sideold1,2.0) + pow(side2-sideold2,2.0));
            fre.row(node0) = fre.row(node0) + k*( (side0-sideold0)*(-vector0/side0) + (side2-sideold2)*(vector2/side2) );
            fre.row(node1) = fre.row(node1) + k*( (side1-sideold1)*(-vector1/side1) + (side0-sideold0)*(vector0/side0) );
            fre.row(node2) = fre.row(node2) + k*( (side2-sideold2)*(-vector2/side2) + (side1-sideold1)*(vector1/side1) );
        }else if ( isDeformArea == true && isInsertionPatch == false ){
            deformnumber_area ++;
            double meanside = sqrt( 4.0*a0/sqrt(3.0) );
            Ere(i) = k/2.0*(pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0));
            fre.row(node0) = fre.row(node0) + k*( (side0-meanside)*(-vector0/side0) + (side2-meanside)*(vector2/side2) );
            fre.row(node1) = fre.row(node1) + k*( (side1-meanside)*(-vector1/side1) + (side0-meanside)*(vector0/side0) );
            fre.row(node2) = fre.row(node2) + k*( (side2-meanside)*(-vector2/side2) + (side1-meanside)*(vector1/side1) );
        }else if ( isDeformShape == true && isDeformArea == false && isInsertionPatch == false ){
            deformnumber_shape ++;
            double meansideold = sqrt( 4.0*area/sqrt(3.0) );
            Ere(i) = k/2.0*(pow(side0-meansideold,2.0) + pow(side1-meansideold,2.0) + pow(side2-meansideold,2.0));
            fre.row(node0) = fre.row(node0) + k*( (side0-meansideold)*(-vector0/side0) + (side2-meansideold)*(vector2/side2) );
            fre.row(node1) = fre.row(node1) + k*( (side1-meansideold)*(-vector1/side1) + (side0-meansideold)*(vector0/side0) );
            fre.row(node2) = fre.row(node2) + k*( (side2-meansideold)*(-vector2/side2) + (side1-meansideold)*(vector1/side1) );
        }else if (isInsertionPatch == true){
            deformnumber_area ++;
            double meanside = sqrt( 4.0*param.s0/sqrt(3.0) );
            //double meanside = param.meanL;
            Ein(i) = K/2.0*(pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0));
            fre.row(node0) = fre.row(node0) + K*( (side0-meanside)*(-vector0/side0) + (side2-meanside)*(vector2/side2) );
            fre.row(node1) = fre.row(node1) + K*( (side1-meanside)*(-vector1/side1) + (side0-meanside)*(vector0/side0) );
            fre.row(node2) = fre.row(node2) + K*( (side2-meanside)*(-vector2/side2) + (side1-meanside)*(vector1/side1) );   
        }
    }
    #pragma omp parallel for reduction(+:E_regularization)
    for (int i = 0; i < face.n_rows; i++){
        E_regularization = E_regularization + Ere(i);
    }
    #pragma omp parallel for reduction(+:E_insertions)
    for (int i = 0; i < face.n_rows; i++){
        E_insertions = E_insertions + Ein(i);
    }
    Fre = fre;
    // note, some sides are considered twice, so as the energy and force is calculated twice, but some sides are not, we just keep the total energy and force since they should vanish finally.
    //E_regularization = E_regularization;
    //Fre = Fre;
    /////////////////////////////////////////////////////////////
    deformnumbers << deformnumber_area << deformnumber_shape << face.n_rows -deformnumber_area-deformnumber_shape << endr;
}

// A: one_ring_vertices; A(i,1)==A(i,2), means face(i) is irregular, has 11 one-ring vertices
void Energy_and_Force(mat vertex, mat vertexold, Mat<int> nearby_vertices, Mat<int> vertexi_face, Mat<int> face, Mat<int> A, Param& param, rowvec spontcurv, Mat<int> insertionpatch, Row<int> isinsertionpatch, Row<double>& E, mat& F, rowvec& deformnumbers, rowvec gqcoeff, cube shape_functions) {
    // initialize the output parameters;
    double numvertex = vertex.n_rows;
    E.fill(0.0);
    F.fill(0.0);
    bool isInsertionAreaConstraint = param.isInsertionAreaConstraint;
    // update the total area and volume;
    rowvec elementS(face.n_rows); elementS.fill(0);
    rowvec elementV(face.n_rows); elementV.fill(0);
    cell_area_volume(vertex,face,A,param.GaussQuadratureN,elementS,elementV,gqcoeff,shape_functions);
    double S = sum(elementS); param.S = S; 
    double V = sum(elementV); param.V = V;
    /////////////////////////////////////////////
    // bending force and energy, constraint force and energy
    //cube fbe(numvertex,3,face.n_rows); fbe.fill(0);  // store bending force and constraint force
    mat fbe(numvertex,3); fbe.fill(0);
    rowvec Ebending(face.n_rows); Ebending.fill(0);
    
    #pragma omp parallel for reduction(+:fbe) 
    for ( int i = 0; i < face.n_rows; i++) {
        // one ring vertices
        mat dots(12,3); 
        for (int j = 0; j < 12; j++) {
            int nodenum = A(i,j);
            dots.row(j) = vertex.row(nodenum);
        }
        ////////////////////////////////////////
        // spontaneous curvature of each patch, 
        double c0 = spontcurv(i);
        ///////////////////////////////////////
        if ( A(i,0) != A(i,1) ) { // for regular patch
            double ebe = 0.0; 
            mat f1(12,3); f1.fill(0); // bending or curvature term
            mat f2(12,3); f2.fill(0); // area term
            mat f3(12,3); f3.fill(0); // vesicle volume term
            element_energy_force_regular(dots, param, c0, ebe, f1, f2, f3, gqcoeff, shape_functions);
            Ebending(i) = ebe;
            for (int j = 0; j < 12; j++) {
                int nodenum = A(i,j);
                fbe.row(nodenum) = fbe.row(nodenum) + f1.row(j) + f2.row(j) + f3.row(j);
            }
        } else {    // A(i,0)==A(i,1) // for irregular patch, one more subdivision is conducted
            double ebe = 0.0; 
            mat f1(11,3); f1.fill(0);// bending or curvature term
            mat f2(11,3); f2.fill(0);// area term
            mat f3(11,3); f3.fill(0);// volume term
            element_energy_force_irregular(i, A, vertex, param, c0, ebe, f1, f2, f3, gqcoeff, shape_functions);
            Ebending(i) = ebe;
            for (int j = 0; j < 11; j++) {
                int nodenum = A(i,j+1);
                fbe.row(nodenum) = fbe.row(nodenum) + f1.row(j) + f2.row(j) + f3.row(j); 
            }
        }
    }
    double E_bending = 0.0;  // curvature 
    #pragma omp parallel for reduction(+:E_bending)
    for ( int i = 0; i < face.n_rows; i++) {
        E_bending = E_bending + Ebending(i);
    }
    mat Fbe = - fbe;     
    
    double us = param.us;
    double uv = param.uv;
    double S0 = param.S0;
    double V0 = param.V0;
    double E_constraint = 0.5*us*pow(S-S0,2.0) + 0.5*uv*pow(V-V0,2.0);
    ///////////////////////////////////////////////////////////////////
    // regularization force and energy
    mat fre(numvertex,3); fre.fill(0.0);
    double E_regularization = 0.0; 
    double E_insertions = 0.0;
    energy_force_regularization(vertex, vertexold, face, param, insertionpatch, isinsertionpatch, E_regularization, fre, E_insertions, deformnumbers);
    ///////////////////////////////////////////////////////////////////////////
    // the total force is the sum of internal, constraint and regularization force.
    mat fxyz  = Fbe + fre;
    F = fxyz;
    ////////////////////////////////////////////////////////////////////////////
    // the total energy is the sum of bending, constraint and regularization energy
    E << E_bending << E_constraint << 0.0 << E_regularization << E_insertions << E_bending + E_constraint + E_regularization + E_insertions << endr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

double line_search_a_to_minimize_energy(double a0, Mat<double>& dx, Mat<double> vertex, Mat<double>& vertexref, Mat<int> vertexi_nearby, Mat<int> vertexi_face, Mat<int> face, Mat<int> face_ring_vertex, Param param, rowvec spontcurv, Mat<int> insertionpatch, Row<int> Isinertionpatch, bool isNCGstucked, rowvec energy0, Mat<double> force0, rowvec gqcoeff, cube shape_functions) {
    double a;
    int numvertex = vertex.n_rows;
    mat force1(numvertex,3);
    rowvec energy1(6);
    rowvec deformnumbers(3);

    a = a0;
    bool isCriterionSatisfied = false;
    double c1 = 1e-4;
    double c2 = 0.1;
    double E0 = energy0(5);
    double E1;
    bool isEnergyDecreased = false;
    bool haveChangedx = false;
    while ( isCriterionSatisfied == false ) {
        mat vertexnew = update_vertex(vertex, a, dx, insertionpatch, face, param);
        Energy_and_Force(vertexnew, vertexref, vertexi_nearby, vertexi_face, face, face_ring_vertex, param, spontcurv, insertionpatch, Isinertionpatch, energy1, force1, deformnumbers, gqcoeff, shape_functions);
        E1 = energy1(5);
        if ( E1 < E0 - 1.0e-20 ) {
            E0 = E1;
            isEnergyDecreased = true;
        }
        if ( isEnergyDecreased == true && E1 > E0 + 1.0e-20 ) {
            isCriterionSatisfied = true;
            a = a/0.8;
            break;
        }
        
        if ( abs(a) < 1e-50 && haveChangedx == false ) {
            cout<<"Note: cannot find an efficient samll a to minimize the energy, now change the direction to gradient!"<<endl;
            dx = force0;
            vertexref = vertex;
            haveChangedx = true;
            a = a0;
        }else if (abs(a) < 1e-50 && haveChangedx == true){
            cout<<"Note: cannot find an efficient samll a to minimize the energy with the direction of gradient!"<<endl;
            a = -1;
            break;
        }
        a = a * 0.8;
    }
    cout<<"value of stepsize: initial a = "<<a0<<", final a = "<<a<<endl;
    return a;
}

void read_struture_vertex(Mat<double>& vertex, char* filename) {
    ifstream fin(filename);
    string line;
    int i = 0;
    while ( getline(fin,line) ) {
        istringstream sin(line);
        vector<string> positions;
        string info;
        while (getline(sin, info, ',')) {
            positions.push_back(info);
        }
        string xstr = positions[0];
        string ystr = positions[1];
        string zstr = positions[2];
        double x, y, z;
        stringstream sx, sy, sz;
        sx << xstr;
        sy << ystr;
        sz << zstr;
        sx >> x;
        sy >> y;
        sz >> z;
        vertex(i,0) = x;
        vertex(i,1) = y;
        vertex(i,2) = z;
        i++;
    }
    if ( i != vertex.n_rows ) {
        cout<< "Wrong! vertices number is "<<vertex.n_rows<<" but read number i = "<< i << endl;
    }
}

void printoutGlobalNCG(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_Global_NCG.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}
void printoutGlobalGD(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_Global_GD.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}
void printoutstuck(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_stuck.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printoutREF(mat vertex1){
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_reference.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printoutforce(mat vertex1){
    int numvertex = vertex1.n_rows;
    ofstream outfile("force.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printout_spontaneouscurvature(rowvec spontcurv){
    ofstream outfile11("spont.csv");
    for (int i = 0; i < spontcurv.n_cols; i++) {
        outfile11 << setprecision(16) << spontcurv(i)+0.0 << '\n';
    }
    outfile11.close(); 
}
