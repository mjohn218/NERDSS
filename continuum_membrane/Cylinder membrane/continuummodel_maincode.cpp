#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include "functions_file1.cpp"

using namespace std;
using namespace arma;

// basic parameters
double Radius  = 500;            // radius of the cylinder, nm
double Length  = 60.0;            // length of the cylinder, nm
double l   = 1.0;                      // triangular side length, nm
double C0  = -0.35;                      // spontaneous curvature of insertion
double c0  = 0.0;                      // spontaneous curvature of membrane
double ds  = 0.0;                      // insertion area
double kc  = 20*4.17;                  // pN.nm  
double us  = 250.0;                  // pN/nm, area stretching modulus; 0.5*us*(ds)^2/s0;
double k   = (1.0e1)*kc;             // coefficient of the regulerization constraint, 
double K   = 2.0*k;                  // spring constant for insertion zones
double gama_shape = 0.2;
double gama_area = 0.2;
bool   isInsertionAreaConstraint = true;
double sigma = 0.0;              // 2*sigma is the lengthscale of decaying spontaneous curvature
bool   isAdditiveScheme = false; // additve scheme for the expansion of spontaneous curvature
int    GaussQuadratureN = 2; 
int    N   = 1e5;                      // total step of iteration
double criterion_force = 1.0e-2;
double criterion_E = 1.0e-3; 
double criterion_S = 1e-5;
double criterion_Er = 1e-5;
bool   isGlobalConstraint = true;        // whether to use Global constraints at the beginning of the simulation
bool   isBoundaryFixed = false;
bool   isBoundaryPeriodic  = true;
bool   isBoundaryFree = false;
bool   isBoundaryForced = false;
double extforce = 0.0;     // pN
/////////////////////////////////////////////////////////////////////////////////
// main code
int main() {
    srand((unsigned)time(NULL)); 
    Param param;
    param.Radius = Radius; param.Length = Length; param.l = l;
    param.isBoundaryFixed = isBoundaryFixed; param.isBoundaryPeriodic = isBoundaryPeriodic; param.isBoundaryFree = isBoundaryFree;
    param.isBoundaryForced = isBoundaryForced; param.extforce = extforce;
    // build the triangular mesh plane. NOTE: ghost vertices and faces are included. All the boundary faces are ghost faces, which should be eliminated when output faces
    Mat<double> vertex = setvertex_Loop_scheme(param); // vertex position
    Mat<int> face = setface_Loop_scheme(param); // face and its surrounding vertex
    Row<int> isEndNode = isEndVertex(param); // element 1 means this vertex is on end
    Row<int> isBoundaryFace = isBoundaryface(face,vertex,param); // element 1 means this face is on boundary
    ///////////////////////////////////////
    // read a structure file
    //char name[32] = "vertex_read.csv";
    //read_struture_vertex(vertex, name);
    ///////////////////////////////////////
    int numvertex = vertex.n_rows;
    int numface = face.n_rows;
    // cout<<numvertex<<", "<<numface<<endl;
    Mat<int> vertexi_face = vertexi_face_with_it(vertex, face); // find out what faces that have vertex_i, 6 or 3 or 2 or 1.
                                                                // if vertexi_face(i, 3-5) == -1, then this vertex i has only 3 faces nearby
    Mat<int> vertexi_nearby = neighbor_vertices(vertexi_face, face); // find out the nearby vertices, 6 or 5(if vertexi_nearby(i,5)==-1)
    Mat<int> face_ring_vertex = one_ring_vertices(face, vertexi_face, isBoundaryFace); // find out the ring_vertices, 12 or 11(if ).
    //////////////////////////////////////////////////////
    // output the vertex and face matrix
    ofstream outfile("face.csv");
    for (int i = 0; i < numface; i++) {
        //if ( isBoundaryFace(i) == 1 ) continue; // ignore the ghost faces.
        outfile << face(i,0)+1 << ',' << face(i,1)+1 << ',' << face(i,2)+1 << '\n';
    }
    outfile.close();
    ofstream outfile1("vertex_begin.csv");
    for (int i = 0; i < numvertex; i++) {
        outfile1 << setprecision(16) << vertex(i,0)+0.0 << ',' << vertex(i,1)+0.0 << ',' << vertex(i,2)+0.0 << '\n';
    }
    outfile1.close();
    /*
    ofstream outfile2("isBoundaryFace.csv");
    for (int i = 0; i < numface; i++) {
        outfile2 << isBoundaryFace(i) << '\n';
    }
    outfile2.close();
    */
    ///////////////////////////////////////////////////////
    Mat<int> insertionpatch; 
             insertionpatch << 105227 << 105228 << 111512 << 111513 << endr; // L=60nm
             insertionpatch = insertionpatch -1;              
    Row<int> Isinsertionpatch = determine_isinsertionpatch(face,insertionpatch);
    int zonenumber = insertionpatch.n_rows; // the number of insertion zones.
    //////////////////////////////////////////////////////////
    // gauss_quadrature and shape functions
    mat VWU = setVMU(GaussQuadratureN);
    cube shape_functions(12,7,VWU.n_rows);
    for (int i = 0; i < VWU.n_rows; i++) {
        rowvec vwu = VWU.row(i);
        mat sf(12,7);
        shapefunctions(vwu,sf);          // 12 shape functions
        shape_functions.slice(i) = sf;
    }
    rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN);
    //////////////////////////////////////////////////////////
    rowvec elementS(face.n_rows); elementS.fill(0);
    membrane_area(vertex,face,isBoundaryFace,face_ring_vertex,GaussQuadratureN,elementS,gqcoeff,shape_functions); // calculate the elemental area 
    double S0 = sum(elementS); // total area
    cout <<"The total area = "<<S0<<endl;
    param.kc = kc; param.us = us/S0; param.k = k; param.K = K; 
    param.C0 = C0; param.c0 = c0; param.gama_shape = gama_shape; param.gama_area = gama_area; param.sigma = sigma; param.GaussQuadratureN = GaussQuadratureN;
    param.isInsertionAreaConstraint =isInsertionAreaConstraint;
    param.isAdditiveScheme = isAdditiveScheme; 
    param.s0 = 2.0/4.0; // /insertionpatch.n_cols; 
    param.S0 = S0; 
    ///////////////////////////////////////
    Row<double> energy(6); energy.fill(0.0); // E_bending, E_constraint, E_memvolume, E_regularization, E_insert, E_tot
    mat force(numvertex,3); force.fill(0);
    mat vertexref = vertex; 
    rowvec deformnumbers(3);
    rowvec spontcurv = determine_spontaneous_curvature(param, face, vertex, insertionpatch);
    printout_spontaneouscurvature(spontcurv);  
    
    //check_nodal_force(vertex, vertexref, vertexi_nearby, vertexi_face, face, isEndNode, isBoundaryFace, face_ring_vertex, param, spontcurv, insertionpatch, Isinsertionpatch, gqcoeff, shape_functions);
    
    Energy_and_Force(vertex, vertexref, vertexi_nearby, vertexi_face, face, isBoundaryFace, face_ring_vertex, param, spontcurv, insertionpatch, Isinsertionpatch, energy, force, deformnumbers, gqcoeff, shape_functions);

    mat Energy(N,6); Energy.fill(0); Energy.row(0) = energy;
    vec forcescale = force_scale(force);
    vec MeanForce(N); MeanForce.fill(0); MeanForce(0) = mean(forcescale);
    vec totalarea(N); totalarea(0) = param.S;

    mat force0 = force; 
    mat force1(force.n_rows,3);
    mat vertex0 = vertex;
    mat vertex1(vertex.n_rows,3);
    mat s0 = force0; 
    mat s1(force.n_rows,3); 
    rowvec energy0 = energy;
    rowvec energy1(6);
    double a0;
    double a1;
    
    bool isNCGstucked = false;
    bool isMinimized = false;
    bool isCriteriaSatisfied = false; 
    bool updateReference = true;
    int i = 0;
    while ( isCriteriaSatisfied == false && i < N-1){
       // updates 
       if ( i == 1 || i%50 == 0 ){
           a0 = l/max(force_scale(s0)); //a0 = 1; 
       }else{
           a0 = a0 * 2e1;
       } 
       {
           a1 = line_search_a_to_minimize_energy(a0, s0, vertex0, vertexref, vertexi_nearby, vertexi_face, face, isBoundaryFace, face_ring_vertex, param, spontcurv, insertionpatch, Isinsertionpatch, isNCGstucked, energy0, force0, gqcoeff, shape_functions);
           vertex1 = update_vertex(vertex0, face, a1, s0, param, isBoundaryFace);
           if ( a1 == -1 ){
               cout<<"step: "<<i<<". Note: no efficent step size a is found. Stop now! "<<endl;
               printoutREF(vertexref);
               printoutstuck(vertex0);
               break;
           }
           // calculate the new force and energy
           Energy_and_Force(vertex1, vertexref, vertexi_nearby, vertexi_face, face, isBoundaryFace, face_ring_vertex, param, spontcurv, insertionpatch, Isinsertionpatch, energy1, force1, deformnumbers, gqcoeff, shape_functions);
           // calculate the direction s
           rowvec df0 = tovector(-force0);
           rowvec df1 = tovector(-force1);
           mat shu0 = df0 * strans(df0);
           mat shu1 = df1 * strans(df1); 
           mat shu10 = df1 * strans(df0);
           double peta1= shu1(0,0) / shu0(0,0); // for local-area constraints, if NCG can't work anylonger, turn to simple gradient method.
           //double peta1 = ( shu1(0,0) - shu10(0,0) ) / shu0(0,0); if ( peta1 < 0 || isNCGstucked == true ){ peta1 = 0; }
           s1 = force1 + peta1 * s0;
           // update the parameters
           vertex0 = vertex1;
           force0 = force1; 
           s0 = s1;
           a0 = a1;
           energy0 = energy1;
           
           i=i+1;
           // store the energy and nodal force
           Energy.row(i) = energy1;
           vec forcescale = force_scale(force1);
           MeanForce(i) = mean(forcescale); 
           totalarea(i) = param.S;
           if (updateReference == true){
               if ( abs(Energy(i,5)-Energy(i-1,5)) < 1e-3 || abs(MeanForce(i)-MeanForce(i-1)) < 1e-3 ){
                   cout<<"update the reference structure! "<<endl;
                   vertexref = vertex1;
               }
           }
       }
       // output parameters     
       cout<<"step: "<< i <<". Sratio= "<<totalarea(i)/S0<<". Energy= "<<Energy(i,5)<<". meanF= "<<MeanForce(i)<<". a= "<<a0<<endl;
       cout<<"step: "<< i <<". Deform number: Area = "<<deformnumbers(0)<<", Shape = "<<deformnumbers(1)<<", nodeform = "<<deformnumbers(2)<<endl;
       // check whether to stop. if the total energy is flat for 100 simulation steps, then stop
       if ( isBoundaryFixed == true || isBoundaryForced == true ){
           if ( i > 500 && abs((Energy(i,5)-Energy(i-500,5))/500) < criterion_E ){
               cout<<"The energy is minimized. Stop now!"<<endl;
               isCriteriaSatisfied = true;
               break;
           }
       }else{
            if ( i > 500 && abs((Energy(i,5)-Energy(i-500,5))/500) < criterion_E && MeanForce(i) < criterion_force ){
               cout<<"The energy is minimized. Stop now!"<<endl;
               isCriteriaSatisfied = true;
               break;
           }
       }
       double Ere_vs_Etot = (Energy(i,3)+0.0)/(Energy(i,5)+0.0);

       ///////////////////////////////////////////////////////
       // output the newvertex
       /*
        if ( i%100 == 0 ){    
            int kk = i/100;
            char filename[20] = "vertex%d.csv";
            sprintf(filename,"vertex%d.csv",kk);
            ofstream outfile(filename);
            for (int j = 0; j < vertex1.n_rows; j++) {
                outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
            }
            outfile.close();
        }
        */
        ofstream outfile2("energy_force.csv"); // E_bending, E_constraint, E_regularization, E_tot, Ere/Etot, s/s0, v/v0, meanForce
        for (int j = 0; j <= i; j++) {
            Ere_vs_Etot = (Energy(j,3)+0.0)/(Energy(j,5)+0.0);
            outfile2 << Energy(j,0)+0.0 << ',' << Energy(j,1)+0.0 << ',' << Energy(j,2)+0.0 << ',' << Energy(j,3)+0.0 << ',' << Energy(j,4)+0.0 << ',' << Energy(j,5)+0.0 << ',' << Ere_vs_Etot << ',' << totalarea(j)/S0 << ',' << ',' << MeanForce(j) << '\n';
        }
        outfile2.close();
    }
    // output the final structure
    ofstream outfile3("vertex_final.csv");
    for (int j = 0; j < vertex0.n_rows; j++) {
        outfile3 << setprecision(16) << vertex0(j,0)+0.0 << ',' << vertex0(j,1)+0.0 << ',' << vertex0(j,2)+0.0 << '\n';
    }
    outfile3.close();
}
