#include <math.h>
#include <iostream>
#include <omp.h>
//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "functions_file2.cpp"

void element_energy_force_regular(vector<vector<double>>& dots, Param& param, double c0, double& Hmean, vector<double>& normVector, double& Ebe, vector<vector<double>>& F_be, vector<vector<double>>& F_s, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions) {
    // F_be is the force related to the curvature; F_s is the force related to the area-constraint; F_v is the force related to the volume-constraint.
    // initialize output parameters
    Ebe = 0.0;
    Hmean = 0.0;
    //////////////////////////////////////////////////////////////
    double kc = param.kc;
    double us = param.us/param.S0;
    double S0 = param.S0;
    double S  = param.S;
    // Gaussian quadrature, second-order or 3 points.
    for (int i = 0; i < GaussQuadratureCoeff.size(); i++) {
        vector<vector<double>> sf = ShapeFunctions[i].sf; //12 shape functions
        // a_1,2,3 covariant vectors; a1,2 contravariant vectors;
        // a_11: a_1 differential to v; a_12: a_1 differential to w;
        vector<double> x = sf[0] * dots;
        vector<double> a_1 = sf[1] * dots;
        vector<double> a_2 = sf[2] * dots;
        vector<double> a_11 = sf[3] * dots;
        vector<double> a_22 = sf[4] * dots;
        vector<double> a_12 = sf[5] * dots;
        vector<double> a_21 = sf[6] * dots;
        vector<double> xa = cross(a_1,a_2);
        double sqa = norm(xa);
        vector<double> xa_1 = cross(a_11,a_2) + cross(a_1,a_21);
        vector<double> xa_2 = cross(a_12,a_2) + cross(a_1,a_22);
        double sqa_1 = 1.0/sqa * dot(xa, xa_1);
        double sqa_2 = 1.0/sqa * dot(xa, xa_2);
        vector<double> a_3 = xa/sqa;
        vector<double> a_31 = 1.0/sqa/sqa * (xa_1*sqa - xa*sqa_1);
        vector<double> a_32 = 1.0/sqa/sqa * (xa_2*sqa - xa*sqa_2);
        vector<double> d = a_3;
        vector<double> d_1 = a_31;
        vector<double> d_2 = a_32;
        vector<double> a1 = cross(a_2,a_3)/sqa;
        vector<double> a2 = cross(a_3,a_1)/sqa;
        vector<double> a11 = 1.0/sqa/sqa *( (cross(a_21,a_3)+cross(a_2,a_31))*sqa - cross(a_2,a_3)*sqa_1 );
        vector<double> a12 = 1.0/sqa/sqa *( (cross(a_22,a_3)+cross(a_2,a_32))*sqa - cross(a_2,a_3)*sqa_2 );
        vector<double> a21 = 1.0/sqa/sqa *( (cross(a_31,a_1)+cross(a_3,a_11))*sqa - cross(a_3,a_1)*sqa_1 );
        vector<double> a22 = 1.0/sqa/sqa *( (cross(a_32,a_1)+cross(a_3,a_12))*sqa - cross(a_3,a_1)*sqa_2 );
        /*
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
        */
        double H_curv = 0.5 * ( dot(a1,d_1) + dot(a2,d_2) );
        vector<double> n1_be = -kc*(2.0*H_curv-c0)*(dot(a1,a1)*d_1 + dot(a1,a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
        vector<double> n2_be = -kc*(2.0*H_curv-c0)*(dot(a2,a1)*d_1 + dot(a2,a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
        vector<double> m1_be = kc*(2.0*H_curv-c0)*a1;
        vector<double> m2_be = kc*(2.0*H_curv-c0)*a2;
        vector<double> n1_cons = us * (S - S0)*a1;
        vector<double> n2_cons = us*(S-S0)*a2;

        double ebe = 0.5 * kc * sqa * pow(2.0*H_curv-c0,2.0);    // bending energy

        vector<vector<double>> f_be(12, vector<double>(3,0.0)); 
        vector<vector<double>> f_cons(12, vector<double>(3,0.0)); 
        for (int j = 0; j < 12; j++) {
            vector<vector<double>> da1 = -sf[3][j]*kron(a1,d) - sf[1][j]*kron(a11,d) - sf[1][j]*kron(a1,d_1) - sf[6][j]*kron(a2,d) - sf[2][j]*kron(a21,d) - sf[2][j]*kron(a2,d_1);
            vector<vector<double>> da2 = -sf[5][j]*kron(a1,d) - sf[1][j]*kron(a12,d) - sf[1][j]*kron(a1,d_2) - sf[4][j]*kron(a2,d) - sf[2][j]*kron(a22,d) - sf[2][j]*kron(a2,d_2);
            vector<double> tempf_be = n1_be*sf[1][j] + m1_be*da1 + n2_be*sf[2][j] + m2_be*da2;
            f_be[j] = - tempf_be * sqa; // the force is the negative derivative 
            vector<double> tempf_cons = n1_cons*sf[1][j] + n2_cons*sf[2][j];
            f_cons[j] = - tempf_cons*sqa; // the force is the negative derivative 
        }
        //Hmean = Hmean + 1.0/2.0*GaussQuadratureCoeff[i]*H_curv;
        //Ebe   = Ebe   + 1.0/2.0*GaussQuadratureCoeff[i]*ebe;
        //F_be  = F_be  + 1.0/2.0*GaussQuadratureCoeff[i]*f_be;
        //F_s   = F_s   + 1.0/2.0*GaussQuadratureCoeff[i]*f_cons;
        Hmean += 1.0/2.0*GaussQuadratureCoeff[i]*H_curv;
        Ebe   += 1.0/2.0*GaussQuadratureCoeff[i]*ebe;
        F_be  += 1.0/2.0*GaussQuadratureCoeff[i]*f_be;
        F_s   += 1.0/2.0*GaussQuadratureCoeff[i]*f_cons;
        normVector += 1.0/2.0*GaussQuadratureCoeff[i]*d;
    }
}

// rgularization energy and forces
void energy_force_regularization(vector<Vertex>& vertex, vector<Face>& face, Param& param){
    double k  = param.k;
    double K  = param.K;
    int deformnumber_shape = 0;
    int deformnumber_area = 0;
    vector<vector<double>> fre(vertex.size(), vector<double>(3,0.0));
    #pragma omp parallel for reduction(+:fre)
    for (int i = 0; i < face.size(); i++){    
        bool IsInsertionPatch = face[i].IsInsertionPatch;
        double Ere = 0.0;

        int node0 = face[i].AdjacentVertex[0]; // three nodes of this face element
        int node1 = face[i].AdjacentVertex[1];
        int node2 = face[i].AdjacentVertex[2];
        vector<double> vector0 = vertex[node0].Coord - vertex[node1].Coord;  double side0 = norm(vector0);
        vector<double> vector1 = vertex[node1].Coord - vertex[node2].Coord;  double side1 = norm(vector1);
        vector<double> vector2 = vertex[node2].Coord - vertex[node0].Coord;  double side2 = norm(vector2);
        double s = (side0 + side1 + side2)/2.0;
        double area = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
        double meanside = 1.0/3.0*(side0 + side1+ side2);
        double gama = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
        vector<double> vectorold0 = vertex[node0].ReferenceCoord - vertex[node1].ReferenceCoord;  double sideold0 = norm(vectorold0);
        vector<double> vectorold1 = vertex[node1].ReferenceCoord - vertex[node2].ReferenceCoord;  double sideold1 = norm(vectorold1);
        vector<double> vectorold2 = vertex[node2].ReferenceCoord - vertex[node0].ReferenceCoord;  double sideold2 = norm(vectorold2);
        s = (sideold0 + sideold1 + sideold2)/2.0;
        double areaold = sqrt( s*(s-sideold0)*(s-sideold1)*(s-sideold2) ); //double areaold = S0/face.n_rows;
            
        bool isDeformShape = false;
        if ( gama > param.gama_shape && param.usingRpi == true ){
            isDeformShape = true;
        }
        bool isDeformArea = false;
        double a0 = areaold; 
        if ( abs(area-a0)/a0 >= param.gama_area && param.usingRpi == true ){
            isDeformArea = true;
        }
            
        if ( isDeformShape == false && isDeformArea == false ){
            Ere = k/2.0*(pow(side0-sideold0,2.0) + pow(side1-sideold1,2.0) + pow(side2-sideold2,2.0));
            fre[node0] += k*( (side0-sideold0)*(-vector0/side0) + (side2-sideold2)*(vector2/side2) );
            fre[node1] += k*( (side1-sideold1)*(-vector1/side1) + (side0-sideold0)*(vector0/side0) );
            fre[node2] += k*( (side2-sideold2)*(-vector2/side2) + (side1-sideold1)*(vector1/side1) );
        }else if ( isDeformArea == true ){
            deformnumber_area ++;
            double meanside = sqrt( 4.0*a0/sqrt(3.0) );
            Ere = k/2.0*(pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0));
            fre[node0] += k*( (side0-meanside)*(-vector0/side0) + (side2-meanside)*(vector2/side2) );
            fre[node1] += k*( (side1-meanside)*(-vector1/side1) + (side0-meanside)*(vector0/side0) );
            fre[node2] += k*( (side2-meanside)*(-vector2/side2) + (side1-meanside)*(vector1/side1) );
        }else if ( isDeformShape == true && isDeformArea == false ){
            deformnumber_shape ++;
            double meansideold = sqrt( 4.0*area/sqrt(3.0) );
            Ere = k/2.0*(pow(side0-meansideold,2.0) + pow(side1-meansideold,2.0) + pow(side2-meansideold,2.0));
            fre[node0] += k*( (side0-meansideold)*(-vector0/side0) + (side2-meansideold)*(vector2/side2) );
            fre[node1] += k*( (side1-meansideold)*(-vector1/side1) + (side0-meansideold)*(vector0/side0) );
            fre[node2] += k*( (side2-meansideold)*(-vector2/side2) + (side1-meansideold)*(vector1/side1) );
        }
        // store energy
        face[i].Energy.EnergyRegularization = Ere;
    }
    ///////////////////////////////////////////////////
    // store force
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].Force.ForceRegularization = fre[i];
    }
    param.DeformationList.DeformShape     = deformnumber_shape;
    param.DeformationList.DeformArea      = deformnumber_area;
    param.DeformationList.Undeform        = face.size() - deformnumber_shape - deformnumber_area;
}

void Energy_and_Force(vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions) {
    // update the total area and volume;
    calculate_element_area_volume(vertex, face, GaussQuadratureCoeff, ShapeFunctions);
    param.S = sum_Membrane_Area(face);    
    double V = sum_Membrane_Volume(face); // total volume   
    /////////////////////////////////////////////
    // bending force and energy, constraint force and energy
    clear_forceONvertex_and_energyONface(vertex, face);
    
    vector<force> Force(vertex.size());
    #pragma omp parallel for reduction(+:Force) 
    for ( int i = 0; i < face.size(); i++) {
        if ( face[i].IsBoundary == true )
            continue;
        
        // one ring vertices
        vector<vector<double>> dots(12,vector<double>(3,0.0)); 
        for (int j = 0; j < 12; j++) {
            int nodenum = face[i].OneRingVertex[j];
            dots[j] = vertex[nodenum].Coord;
        }
        ////////////////////////////////////////
        // spontaneous curvature of each patch, 
        double c0 = face[i].SpontCurvature;
        { // for regular patch
            double ebe = 0.0;    // curvature energy of this element;
            double h_mean = 0.0; // mean curvature of this element;
            vector<double> normVector(3,0.0);
            vector<vector<double>> F_be(12,vector<double>(3,0.0)); // bending or curvature term
            vector<vector<double>> F_consA(12,vector<double>(3,0.0)); // area term
            element_energy_force_regular(dots, param, c0, h_mean, normVector, ebe, F_be, F_consA, GaussQuadratureCoeff, ShapeFunctions);
            face[i].Energy.EnergyCurvature = ebe;
            face[i].MeanCurvature = h_mean;
            face[i].normVector = normVector;
            //cout<<i<<", "<<h_mean<<", "<<ebe<<endl;    
            for (int j = 0; j < 12; j++) {
                int nodenum = face[i].OneRingVertex[j];
                //vertex[nodenum].Force.ForceCurvature += F_be[j];
                //vertex[nodenum].Force.ForceArea      += F_consA[j];
                force Forcetmp; Forcetmp.ForceCurvature = F_be[j];
                                Forcetmp.ForceArea      = F_consA[j];
                Force[nodenum] += Forcetmp;
            }
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].Force = Force[i];
        vertex[i].Force.ForceTotal = vertex[i].Force.ForceCurvature + vertex[i].Force.ForceArea + vertex[i].Force.ForceVolume + vertex[i].Force.ForceThickness + vertex[i].Force.ForceTilt + vertex[i].Force.ForceRegularization;
    }
    ///////////////////////////////////////////////////////////////////
    // regularization force and energy
    energy_force_regularization(vertex, face, param);
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    // On each vertex, the total force is the sum of internal, constraint and regularization force.
    //#pragma omp parallel for
    //for (int i = 0; i < vertex.size(); i++){
        //vertex[i].Force.ForceTotal = vertex[i].Force.ForceCurvature + vertex[i].Force.ForceArea + vertex[i].Force.ForceVolume + vertex[i].Force.ForceThickness + vertex[i].Force.ForceTilt + vertex[i].Force.ForceRegularization;
    //}
    // On each face, the total energy is the sum of internal, constraint and regularization energy.
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        face[i].Energy.EnergyTotal = face[i].Energy.EnergyCurvature + face[i].Energy.EnergyArea + face[i].Energy.EnergyVolume + face[i].Energy.EnergyThickness + face[i].Energy.EnergyTilt + face[i].Energy.EnergyRegularization;
    }
    //////////////////////////////////////////////////////////////////
    // total energy to store in param.Energy
    //calculate_total_energy_save_in_paramEnergy(face, param);
    energy tmpEnergy;
    #pragma omp parallel for reduction(+:tmpEnergy)
    for (int i = 0; i < face.size(); i++){
        /*
        tmpEnergy.EnergyCurvature      += face[i].Energy.EnergyCurvature;
        tmpEnergy.EnergyArea           += face[i].Energy.EnergyArea;
        tmpEnergy.EnergyVolume         += face[i].Energy.EnergyVolume;
        tmpEnergy.EnergyThickness      += face[i].Energy.EnergyThickness;
        tmpEnergy.EnergyTilt           += face[i].Energy.EnergyTilt;
        tmpEnergy.EnergyRegularization += face[i].Energy.EnergyRegularization;
        */
        tmpEnergy += face[i].Energy;
    }
    double us = param.us; double S0 = param.S0; double S  = param.S;
    if (param.isGlobalConstraint == true) 
        tmpEnergy.EnergyArea = 0.5*us/S0*pow(S-S0,2.0);
    tmpEnergy.EnergyTotal = tmpEnergy.EnergyCurvature + tmpEnergy.EnergyArea + tmpEnergy.EnergyVolume + tmpEnergy.EnergyThickness + tmpEnergy.EnergyTilt + tmpEnergy.EnergyRegularization;
    param.Energy = tmpEnergy;
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    manage_force_for_boundary_ghost_vertex(vertex, face, param);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
double LineSearch_for_StepSize_to_minimize_energy(double a0, vector<vector<double>>& s0, vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions){
    double a = a0;
    int    numvertex = vertex.size();
    double E0 = param.Energy.EnergyTotal;
    double E1;
    //param.isNCGstucked = false;
    param.usingRpi = true;
    
    double NCGfactor0 = 0.0;
    #pragma omp parallel for reduction(+:NCGfactor0) 
    for (int i = 0; i < vertex.size(); i++){
        NCGfactor0 += dot(-vertex[i].ForcePrevious.ForceTotal, s0[i]);
    }
    double c1 = 1e-4; // two factors for the nonlinear conjugate gradient method
    double c2 = 0.1;
    bool   isCriterionSatisfied = false;
    while ( isCriterionSatisfied == false ) {
        a = a * 0.8;
        update_vertex_from_NonlinearConjugateGradient_s0(vertex, face, a, s0, param); // use CoordPrevious to calculate the new coord
        Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
        E1 = param.Energy.EnergyTotal;
        
        if ( param.usingNCG == true && param.isNCGstucked == false ){
            double NCGfactor = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor) 
            for (int i = 0; i < vertex.size(); i++){
                NCGfactor += dot(-vertex[i].Force.ForceTotal, s0[i]);
            }
            //if ( E1 <= E0 + c1*a*NCGfactor0 && NCGfactor >= c2*NCGfactor0 ){ // Wolfe conditions
            if ( E1 <= E0 + c1*a*NCGfactor0 && abs(NCGfactor) <= c2*abs(NCGfactor0) ){ // strong Wolfe conditions
                break;
            }
            if ( a < 1.0e-9 ) {
                cout<<"Now change the NCG WolfeConditions to simple line search method!"<<endl;
                //WolfeConditions = false;
                a = a0;
                update_reference_from_Coord(vertex);
                //isNCGstucked = true;
                param.isNCGstucked = true;
                param.usingRpi = false;
                // change s0 to nodalForce
                #pragma omp parallel for 
                for (int i = 0; i < vertex.size(); i++){
                    s0[i] = vertex[i].ForcePrevious.ForceTotal;
                }
            }
        }else{
            param.isNCGstucked = true;
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

/*
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
*/
/*
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
*/

void check_nodal_force(vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions){
    cout<<"check if the nodal force is correct: "<<endl;
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
    double E = param.Energy.EnergyTotal;
    int numvertex = vertex.size();
    vector<vector<double>> ForceReal(numvertex,vector<double>(3,0.0)); 
    double dx = 1.0e-8;
    for (int i = 0; i < numvertex; i++) {
        if (vertex[i].IsGhost == true) 
            continue;

        vector<Vertex> vertextry = vertex;
        vertextry[i].Coord[0] = vertex[i].Coord[0] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions);
        double Etmpx = param.Energy.EnergyTotal;
        double fx = - (Etmpx - E)/dx;

        vertextry = vertex;
        vertextry[i].Coord[1] = vertex[i].Coord[1] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions);
        double Etmpy = param.Energy.EnergyTotal;
        double fy = - (Etmpy - E)/dx;

        vertextry = vertex;
        vertextry[i].Coord[2] = vertex[i].Coord[2] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions);
        double Etmpz = param.Energy.EnergyTotal;
        double fz = - (Etmpz - E)/dx;

        ForceReal[i][0] = fx;
        ForceReal[i][1] = fy;
        ForceReal[i][2] = fz;
        double Freal = norm(ForceReal[i]);
        double F = norm(vertex[i].Force.ForceTotal);
        double dF = norm(ForceReal[i] - vertex[i].Force.ForceTotal);
        cout<<i<<setprecision(15)<<". Freal = "<<Freal<<", F = "<<F<<", dF = "<<dF<<endl;
    }
}
