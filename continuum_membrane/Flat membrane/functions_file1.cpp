#include <math.h>
#include <iostream>
#include <armadillo>
#include <omp.h>
//#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#include "functions_file2.cpp"

void element_energy_force_regular(mat dots, Param param, double c0, double& Ebe, mat& F_be, mat& F_s, rowvec gqcoeff, cube shape_functions) {
    // F_be is the force related to the curvature; F_s is the force related to the area-constraint; F_v is the force related to the volume-constraint.
    // initialize output parameters
    Ebe = 0.0;
    F_be.fill(0);
    F_s.fill(0);
    //////////////////////////////////////////////////////////////
    double kc = param.kc;
    double us = param.us/param.S0;
    double S0 = param.S0;
    double S  = param.S;
    int GaussQuadratureN = param.GaussQuadratureN;
    // Gaussian quadrature, second-order or 3 points.
    mat VWU = setVMU(GaussQuadratureN); 
    //rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    for (int i = 0; i < VWU.n_rows; i++) {
        double ebe = 0.0;
        mat f_be(12,3); f_be.fill(0);
        mat f_cons(12,3); f_cons.fill(0);
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
        
        rowvec zaxis(3); zaxis << 0.0 << 0.0 << 1.0 << endr;
        mat shu = d * strans(zaxis);
        bool isCurvPositive = false;
        // the triangular norm is towards up. Then the curvature is defined as positive
        if ( shu(0,0) > 0 ){ 
             isCurvPositive = true;
        }
        shu = a1*strans(d_1) + a2*strans(d_2);
        double H_curv;
        rowvec n1_be(3); rowvec n2_be(3); rowvec m1_be(3); rowvec m2_be(3);
        if (isCurvPositive == true){ // For upword patch, the curvature is positive
            H_curv = 0.5*shu(0,0);
            n1_be = -kc*(2.0*H_curv-c0)*(a1*strans(a1)*d_1+a1*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
            n2_be = -kc*(2.0*H_curv-c0)*(a2*strans(a1)*d_1+a2*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
            m1_be = kc*(2.0*H_curv-c0)*a1;
            m2_be = kc*(2.0*H_curv-c0)*a2;
        }else{ // For downword patch, the curvature is negative
            H_curv = - 0.5*shu(0,0);
            n1_be = kc*(2.0*H_curv-c0)*(a1*strans(a1)*d_1+a1*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
            n2_be = kc*(2.0*H_curv-c0)*(a2*strans(a1)*d_1+a2*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
            m1_be = -kc*(2.0*H_curv-c0)*a1;
            m2_be = -kc*(2.0*H_curv-c0)*a2;
        }
        ebe = 0.5*kc*sqa*pow(2.0*H_curv-c0,2.0);    // bending energy
        rowvec n1_cons = us*(S-S0)*a1;
        rowvec n2_cons = us*(S-S0)*a2;
        for (int j = 0; j < 12; j++) {
            mat da1 = -sf(j,3)*kron(strans(a1),d) - sf(j,1)*kron(strans(a11),d) - sf(j,1)*kron(strans(a1),d_1) - sf(j,6)*kron(strans(a2),d) - sf(j,2)*kron(strans(a21),d) - sf(j,2)*kron(strans(a2),d_1);
            mat da2 = -sf(j,5)*kron(strans(a1),d) - sf(j,1)*kron(strans(a12),d) - sf(j,1)*kron(strans(a1),d_2) - sf(j,4)*kron(strans(a2),d) - sf(j,2)*kron(strans(a22),d) - sf(j,2)*kron(strans(a2),d_2);
            rowvec tempf_be = n1_be*sf(j,1) + m1_be*da1 + n2_be*sf(j,2) + m2_be*da2;
            f_be.row(j) = tempf_be*sqa;
            rowvec tempf_cons = n1_cons*sf(j,1) + n2_cons*sf(j,2);
            f_cons.row(j) = tempf_cons*sqa;
        }
        Ebe = Ebe + 1.0/2.0*gqcoeff(i)*ebe;
        F_be = F_be + 1.0/2.0*gqcoeff(i)*f_be;
        F_s = F_s + 1.0/2.0*gqcoeff(i)*f_cons;
    }
}

void attachpatch_energy_force_regular(mat dots, Param param, double s, double c0, double& Ebe, mat& F_be, mat& F_s, rowvec gqcoeff, cube shape_functions) {
    // F_be is the force related to the curvature; F_s is the force related to the area-constraint; F_v is the force related to the volume-constraint.
    // initialize output parameters
    Ebe = 0.0;
    F_be.fill(0);
    F_s.fill(0);
    //////////////////////////////////////////////////////////////
    double kc = param.kc;
    double s0 = param.s0;
    double us = param.us/s0;
    int GaussQuadratureN = param.GaussQuadratureN;
    // Gaussian quadrature, second-order or 3 points.
    mat VWU = setVMU(GaussQuadratureN); 
    //rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    for (int i = 0; i < VWU.n_rows; i++) {
        double ebe = 0.0;
        mat f_be(12,3); f_be.fill(0);
        mat f_cons(12,3); f_cons.fill(0);
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
        
        rowvec zaxis(3); zaxis << 0.0 << 0.0 << 1.0 << endr;
        mat shu = d * strans(zaxis);
        bool isCurvPositive = false;
        // the triangular norm is towards up. Then the curvature is defined as positive
        if ( shu(0,0) > 0 ){ 
             isCurvPositive = true;
        }
        shu = a1*strans(d_1) + a2*strans(d_2);
        double H_curv;
        rowvec n1_be(3); rowvec n2_be(3); rowvec m1_be(3); rowvec m2_be(3);
        if (isCurvPositive == true){ // For upword patch, the curvature is positive
            H_curv = 0.5*shu(0,0);
            n1_be = -kc*(2.0*H_curv-c0)*(a1*strans(a1)*d_1+a1*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
            n2_be = -kc*(2.0*H_curv-c0)*(a2*strans(a1)*d_1+a2*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
            m1_be = kc*(2.0*H_curv-c0)*a1;
            m2_be = kc*(2.0*H_curv-c0)*a2;
        }else{ // For downword patch, the curvature is negative
            H_curv = - 0.5*shu(0,0);
            n1_be = kc*(2.0*H_curv-c0)*(a1*strans(a1)*d_1+a1*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
            n2_be = kc*(2.0*H_curv-c0)*(a2*strans(a1)*d_1+a2*strans(a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
            m1_be = -kc*(2.0*H_curv-c0)*a1;
            m2_be = -kc*(2.0*H_curv-c0)*a2;
        }
        ebe = 0.5*kc*sqa*pow(2.0*H_curv-c0,2.0);    // bending energy
        rowvec n1_cons = us*(s-s0)*a1;
        rowvec n2_cons = us*(s-s0)*a2;
        for (int j = 0; j < 12; j++) {
            mat da1 = -sf(j,3)*kron(strans(a1),d) - sf(j,1)*kron(strans(a11),d) - sf(j,1)*kron(strans(a1),d_1) - sf(j,6)*kron(strans(a2),d) - sf(j,2)*kron(strans(a21),d) - sf(j,2)*kron(strans(a2),d_1);
            mat da2 = -sf(j,5)*kron(strans(a1),d) - sf(j,1)*kron(strans(a12),d) - sf(j,1)*kron(strans(a1),d_2) - sf(j,4)*kron(strans(a2),d) - sf(j,2)*kron(strans(a22),d) - sf(j,2)*kron(strans(a2),d_2);
            rowvec tempf_be = n1_be*sf(j,1) + m1_be*da1 + n2_be*sf(j,2) + m2_be*da2;
            f_be.row(j) = tempf_be*sqa;
            rowvec tempf_cons = n1_cons*sf(j,1) + n2_cons*sf(j,2);
            f_cons.row(j) = tempf_cons*sqa;
        }
        Ebe = Ebe + 1.0/2.0*gqcoeff(i)*ebe;
        F_be = F_be + 1.0/2.0*gqcoeff(i)*f_be;
        F_s = F_s + 1.0/2.0*gqcoeff(i)*f_cons;
    }
}

void energy_force_regularization(mat vertex, mat vertexold, Mat<int> face, Row<int> isGhostFace, Param param, Row<int> isinsertionpatch, double& E_regularization, mat& Fre, double& E_insertions, rowvec& deformnumbers){
    double k  = param.k;
    double K  = param.K;
    //Fre.fill(0.0);
    E_regularization = 0.0;
    E_insertions = 0.0;
    deformnumbers.fill(0);
    int deformnumber_shape = 0;
    int deformnumber_area = 0;
    rowvec Ere(face.n_rows); Ere.fill(0.0);
    rowvec Ein(face.n_rows); Ein.fill(0.0);
    mat fre(vertex.n_rows,3); fre.fill(0.0);
    #pragma omp parallel for reduction(+:fre)
    for (int i = 0; i < face.n_rows; i++){  
        if ( isGhostFace(i) == 1 ){
            continue;
        }
        bool isInsertionPatch = false;
        if ( isinsertionpatch(i) == 1 && abs(param.C0-param.c0)>1e-9){
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
        if ( gama > param.gama_shape && param.usingRpi == true){
            isDeformShape = true;
        }
        bool isDeformArea = false;
        double a0 = areaold;
        if ( abs(area-a0)/a0 >= param.gama_area && param.usingRpi == true){
            isDeformArea = true;
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

void bind_to_curvPoint(mat dots, rowvec point, Param param, double& E, mat& Force){
    E = 0.0;
    Force.fill(0);
    //////////////////////////////////////////////////////////////
    double k = param.bindCoefficient;
    double r0 = param.bindLength;

    int N=6; double w=3.0/8.0/N; 
    rowvec coeff(7); coeff.fill(w); coeff(0) = 1.0-w*N;
    rowvec x = coeff * dots;
    rowvec rvec = point - x;
    double r = norm(rvec,2.0);
    E = 0.5*k*pow(r-r0,2.0);
    for (int i = 0; i < 7; i++) {
        Force.row(i) = - coeff(i) * k*(r-r0)/r * rvec;
    }
}

// A: one_ring_vertices; A(i,1)==A(i,2), means face(i) is irregular, has 11 one-ring vertices
void Energy_and_Force(mat vertex, mat vertexold, Mat<int> face, Row<int> isBoundaryFace, Row<int> isBoundaryVertex, Row<int> isGhostFace, Row<int> isGhostVertex, 
                      Mat<int> A, Mat<int> closestVertex, Param& param, rowvec spontcurv, Row<int> isinsertionpatch, Row<double>& E, mat& F, rowvec& deformnumbers, rowvec gqcoeff, cube shape_functions) {
    // initialize the output parameters;
    double numvertex = vertex.n_rows;
    E.fill(0.0);
    F.fill(0.0);
    bool isInsertionAreaConstraint = param.isInsertionAreaConstraint;
    // update the total area and volume;
    rowvec elementS(face.n_rows); elementS.fill(0);
    membrane_area(vertex,face,isBoundaryFace,A,param.GaussQuadratureN,elementS,gqcoeff,shape_functions);
    double S = sum(elementS); param.S = S; // area exlcude the attachpatch.
    /////////////////////////////////////////////
    // bending force and energy, constraint force and energy
    mat fbe(numvertex,3); fbe.fill(0.0);
    rowvec Ebending(face.n_rows); Ebending.fill(0);
    
    #pragma omp parallel for reduction(+:fbe) 
    for ( int i = 0; i < face.n_rows; i++) {
        if ( isBoundaryFace(i) == 1 )
            continue;
        
        // one ring vertices
        mat dots(12,3); 
        for (int j = 0; j < 12; j++) {
            int nodenum = A(i,j);
            dots.row(j) = vertex.row(nodenum);
        }
        ////////////////////////////////////////
        // spontaneous curvature of each patch, 
        double c0 = spontcurv(i);
        { // for regular patch
            double ebe = 0.0; 
            mat f1(12,3); f1.fill(0); // bending or curvature term
            mat f2(12,3); f2.fill(0); // area term
            element_energy_force_regular(dots, param, c0, ebe, f1, f2, gqcoeff, shape_functions);
            Ebending(i) = ebe;
            for (int j = 0; j < 12; j++) {
                int nodenum = A(i,j);
                fbe.row(nodenum) = fbe.row(nodenum) + f1.row(j) + f2.row(j);
            }
        }
    }
    double E_bending = 0.0;  // curvature 
    #pragma omp parallel for reduction(+:E_bending)
    for ( int i = 0; i < face.n_rows; i++) {
        E_bending = E_bending + Ebending(i);
    }
    double us = param.us;
    double S0 = param.S0;
    double E_constraint = 0.5*us/S0*pow(S-S0,2.0); 
    ///////////////////////////////////////////////////////////////////
    // regularization force and energy
    mat fre(numvertex,3); fre.fill(0.0);
    double E_regularization = 0.0; 
    double E_insertions = 0.0;
    energy_force_regularization(vertex, vertexold, face, isGhostFace, param, isinsertionpatch, E_regularization, fre, E_insertions, deformnumbers);
    ///////////////////////////////////////////////////////////////////////////
    // the total force is the sum of internal, constraint and regularization force.
    mat fxyz = - fbe;
    fxyz  = manage_ghost_force(fxyz, face, isGhostFace, isGhostVertex, param);
    F = fxyz + fre;
    ////////////////////////////////////////////////////////////////////////////
    // the total energy is the sum of bending, constraint and regularization energy
    E << E_bending << E_constraint << 0 << E_regularization << E_insertions << E_bending + E_constraint + E_regularization + E_insertions << endr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

double line_search_a_to_minimize_energy(double a0, Mat<double>& dx, Mat<double> vertex, Mat<double>& vertexref, Mat<int> face, 
                                        Row<int> isBoundaryFace, Row<int> isBoundaryVertex, Row<int> isGhostFace, Row<int> isGhostVertex, Mat<int> face_ring_vertex, Mat<int> closestVertex,
                                        Param& param, rowvec spontcurv, Row<int> Isinertionpatch, rowvec energy0, Mat<double> force0, rowvec gqcoeff, cube shape_functions) {
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
    bool usingNCG = true;
    bool isNCGstucked = false;
    param.usingRpi = true;
    mat shu0 = tovector(-force0) * strans(tovector(dx));

    while ( isCriterionSatisfied == false ) {
        a = a * 0.8;
        mat vertexnew = update_vertex(vertex, face, a, dx, param, isGhostFace, isGhostVertex);
        Energy_and_Force(vertexnew, vertexref, face, isBoundaryFace, isBoundaryVertex, isGhostFace, isGhostVertex, 
                         face_ring_vertex, closestVertex, param, spontcurv, Isinertionpatch, energy1, force1, deformnumbers, gqcoeff, shape_functions);
        E1 = energy1(5);
        
        if ( param.usingNCG == true && isNCGstucked == false ){
            mat shu1 = tovector(-force1) * strans(tovector(dx));
            //if ( energy1(3) <= energy0(3) + c1*a*shu0(0,0) && shu1(0,0) >= c2*shu0(0,0) ){ // Wolfe conditions
            if ( E1 <= E0 + c1*a*shu0(0,0) && abs(shu1(0,0)) <= c2*abs(shu0(0,0)) ){ // strong Wolfe conditions
                //cout<<"NCG method, a = "<<a<<", Ebe = "<<energy1(0)<<", Econs = "<<energy1(1)<<", Ebar = "<<energy1(2)<<", Ere = "<<energy1(3)<<", Etot = "<<energy1(5)<<endl;
                break;
            }
            if ( a < 1.0e-9 ) {
                cout<<"Now change the NCG WolfeConditions to simple line search method!"<<endl;
                //WolfeConditions = false;
                a = a0;
                vertexref = vertex;
                isNCGstucked = true;
                param.usingRpi = false;
                dx = force0;
            }
        }else{
            if ( E1 < E0 ) {
                //cout<<"Simple method, a = "<<a<<", Ebe = "<<energy1(0)<<", Econs = "<<energy1(1)<<", Ebar = "<<energy1(2)<<", Ere = "<<energy1(3)<<", Etot = "<<energy1(5)<<endl;
                break;
            }
            if ( a < 1.0e-20 ) {
                 cout<<"Note: cannot find an efficient samll a to minimize the energy even with the simple method!"<<endl;
                 return -1;
            }
        }
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

void check_nodal_force(mat vertex, mat vertexref, Mat<int> face, Row<int> isBoundaryFace, Row<int> isBoundaryVertex, Row<int> isGhostFace, Row<int> isGhostVertex,
                       Mat<int> face_ring_vertex, Mat<int> closestVertex, Param& param, rowvec spontcurv, Row<int> Isinsertionpatch, rowvec gqcoeff, cube shape_functions){
    cout<<"check if the nodal force is correct: "<<endl;
    rowvec energy(6);
    energy.fill(0);
    mat force(vertex.n_rows,3);
    force.fill(0);
    rowvec deformnumbers(3);
    Energy_and_Force(vertex, vertexref, face, isBoundaryFace, isBoundaryVertex, isGhostFace, isGhostVertex, face_ring_vertex, closestVertex, param, 
                     spontcurv, Isinsertionpatch, energy, force, deformnumbers, gqcoeff, shape_functions);
    double E = energy(5);
    vec forcescale = force_scale(force);
    int numvertex = vertex.n_rows;
    mat forcereal(numvertex,3); forcereal.fill(0);
    double dx = 1.0e-8;
    mat forceha(vertex.n_rows,3);
    for (int i = 0; i < numvertex; i++) {
        if (isGhostVertex(i) == 1) 
            continue;
        mat vertextry = vertex;
        vertextry(i,0) = vertex(i,0)+dx;
        mat forcetemp(numvertex,3);
        forcetemp.fill(0);
        Energy_and_Force(vertextry, vertexref, face, isBoundaryFace, isBoundaryVertex, isGhostFace, isGhostVertex, face_ring_vertex, closestVertex, param, spontcurv, Isinsertionpatch, energy, forcetemp, deformnumbers, gqcoeff, shape_functions);
        double fx = - (energy(5)-E)/dx;

        vertextry = vertex;
        vertextry(i,1) = vertex(i,1)+dx;
        Energy_and_Force(vertextry, vertexref, face, isBoundaryFace, isBoundaryVertex, isGhostFace, isGhostVertex, face_ring_vertex, closestVertex, param, spontcurv, Isinsertionpatch, energy, forcetemp, deformnumbers, gqcoeff, shape_functions);
        double fy = - (energy(5)-E)/dx;

        vertextry = vertex;
        vertextry(i,2) = vertex(i,2)+dx;
        Energy_and_Force(vertextry, vertexref, face, isBoundaryFace, isBoundaryVertex, isGhostFace, isGhostVertex, face_ring_vertex, closestVertex, param, spontcurv, Isinsertionpatch, energy, forcetemp, deformnumbers, gqcoeff, shape_functions);
        double fz = - (energy(5)-E)/dx;

        forcereal(i,0) = fx;
        forcereal(i,1) = fy;
        forcereal(i,2) = fz;
        double frea = norm(forcereal.row(i),2.0);
        double f = norm(force.row(i),2.0);
        double cha = norm(forcereal.row(i)-force.row(i),2.0);
        cout<<i<<setprecision(15)<<". f= "<<f<<", freal= "<<frea<<", cha= "<<cha<<endl;
    }
    mat dforce = forcereal - force;
    vec dforcescale = force_scale(dforce);
    vec dfratio(numvertex);
    dfratio.fill(0);
    for(int i=0; i<numvertex; i++) {
        if (isGhostVertex(i)==1) continue;
        dfratio(i)=dforcescale(i)/forcescale(i);
    }
    cout<<"Minratio = "<<min(dfratio)<<", Maxratio = "<<max(dfratio)<<endl;
    double small=1; double large=0;
    int smallindex=0; int largeindex = 0;
    for (int i=0; i<numvertex; i++){
        if (isGhostVertex(i)==1) continue;
        if (dfratio(i) < small){
            small = dfratio(i);
            smallindex = i;
        }
        if (dfratio(i) > large){
            large = dfratio(i);
            largeindex = i;
        }
    }
    cout<<"Minindex = "<<smallindex<<", Maxindex = "<<largeindex<<endl;
}
