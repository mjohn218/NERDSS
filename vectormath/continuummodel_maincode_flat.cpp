#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include <numeric> // for accumulate vector
//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <algorithm>    // std::min_element, std::max_element
#include "functions_file1.cpp"

using namespace std;
// basic parameters
double sideX  = 50.0;            // rectangle sidelength x, nm
double sideY  = sideX;            // rectangle sidelength y, nm
double l   = 2.0;                 // triangular side length, nm
double C0  = 1.0;                 // spontaneous curvature of insertion. Towards up is defined as positive
double c0  = 0.0;                      // spontaneous curvature of membrane
double ds  = 0.0;                      // insertion area
double miu = 1.0;                  // area constraint coefficient
double kc  = 20*4.17;                  // pN.nm  
double us  = 250.0;                  // pN/nm, area stretching modulus; 0.5*us*(ds)^2/s0;
double k   = (1.0e1)*kc;             // coefficient of the regulerization constraint, 
double K   = 1.0*k;                  // spring constant for insertion zones
double gama_shape = 0.2;
double gama_area = 0.2;
bool   isInsertionAreaConstraint = true;
double sigma = 0.0;              // 2*sigma is the lengthscale of decaying spontaneous curvature
bool   isAdditiveScheme = false; // additve scheme for the expansion of spontaneous curvature
int    GaussQuadratureN = 2; 
int    N   = 1e5;                      // total step of iteration
double CriterionForForce = 1.0e-2;
double CriterionForEnergy = 1.0e-3; 
double CriterionForArea = 1e-5;
double CriterionForRegularization = 1e-5;
bool   isGlobalConstraint = true;        // whether to use Global constraints at the beginning of the simulation
bool   isBoundaryFixed = false;
bool   isBoundaryPeriodic  = true;
bool   isBoundaryFree = false;
/////////////////////////////////////////////////////////////////////////////////
// main code
int main() {
    srand((unsigned)time(NULL)); 

    // build the triangular mesh plane. NOTE: ghost vertices and faces are included. 
    // All the boundary faces are ghost faces, which should be eliminated when output faces
    vector<Vertex> vertex = setVertex_Loop_scheme(sideX, sideY, l); // vertex position
    vector<Face> face = setFace_Loop_scheme(sideX, sideY, l); // face and its surrounding vertex
    determine_Boundary_vertex_face(sideX, sideY, l, vertex, face); // flag whether this vertex or face is on boundary
    determine_Ghost_vertex_face(sideX, sideY, l, isBoundaryFixed, isBoundaryPeriodic, isBoundaryFree, vertex, face); // flag whether this vertex or face is ghost
    determine_AdjacentFace_for_vertex(vertex, face); // find out what faces that have vertex_i, 6 or 3 or 2 or 1.
    determine_AdjacentVertex_for_vertex(vertex, face); // find out the nearby vertices for each vertex, 6 or less.
    determine_OneRingVertex_for_face(vertex, face); // find out the one-ring-vertices, 12 for flat surface with only regular patch.

    //////////////////////////////////////////////////////
    ///////////////////////////////////////
    // read a structure file
    //char name[32] = "vertex_read.csv";
    //read_struture_vertex(vertex, name);
    ///////////////////////////////////////
    int numvertex = vertex.size();
    int numface = face.size();
    // cout<<numvertex<<", "<<numface<<endl;
    // output the vertex and face matrix
    ofstream outfile("face.csv");
    for (int i = 0; i < numface; i++) {
        outfile << face[i].AdjacentVertex[0] + 1 << ',' << face[i].AdjacentVertex[1] + 1 << ',' << face[i].AdjacentVertex[2] + 1 << '\n';
    }
    outfile.close();
    ofstream outfile1("vertex_begin.csv");
    for (int i = 0; i < numvertex; i++) {
        outfile1 << setprecision(16) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << '\n';
    }
    outfile1.close();
    ///////////////////////////////////////////////////////
    vector<vector<int> > InsertionPatch { {650, 651, 652, 653},
                                          {810, 811, 812, 813} };            
    determine_IsInsertionPatch_for_face(face, InsertionPatch);
    determine_SpontaneousCurvature_for_face(C0, c0, face);
    //////////////////////////////////////////////////////////
    // gauss_quadrature and shape functions
    vector<vector<double>> VWU = setVMU(GaussQuadratureN);
    vector<double> GaussQuadratureCoeff = setVMUcoefficient(GaussQuadratureN);
    vector<Shapefunctions> ShapeFunctions(VWU.size());
    for (int i = 0; i < VWU.size(); i++) {
        vector<double> vwu = VWU[i];
        ShapeFunctions[i].sf = determine_ShapeFunctions(vwu);
    }
    
    //////////////////////////////////////////////////////////
    calculate_element_area_volume(vertex, face, GaussQuadratureCoeff, ShapeFunctions); // calculate the elemental area 
    double S0 = sum_Membrane_Area(face); // total area
    double V0 = sum_Membrane_Volume(face); // total volume

    Param param;
    param.kc = kc; param.us = us; param.k = k; param.K = K; 
    param.C0 = C0; param.c0 = c0; param.gama_shape = gama_shape; param.gama_area = gama_area; param.sigma = sigma; param.GaussQuadratureN = GaussQuadratureN;
    param.isInsertionAreaConstraint =isInsertionAreaConstraint;
    param.isAdditiveScheme = isAdditiveScheme; param.isGlobalConstraint = isGlobalConstraint;
    param.s0 = sqrt(3.0)/4.0*l*l;//2.0/4.0; // /insertionpatch.n_cols; 
    param.S0 = S0*miu; 
    param.sideX = sideX; param.sideY = sideY; param.l = l;
    param.isBoundaryFixed = isBoundaryFixed; 
    param.isBoundaryPeriodic = isBoundaryPeriodic; 
    param.isBoundaryFree = isBoundaryFree;
    param.usingNCG = true;
    param.InsertionPatch = InsertionPatch;
    ///////////////////////////////////////
    update_reference_from_Coord(vertex);
    // check whether the code is correct, especially whether the force is correct!
    //check_nodal_force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
    update_PreviousCoord_for_vertex(vertex);
    update_PreviousForce_for_vertex(vertex);

single_particles_diffusing_on_surface(vertex,face, param);
exit(0);

    vector<double> MeanForce(N,0.0); 
    vector<double> AreaTotal(N,0.0); AreaTotal[0] = param.S;
    vector<energy> Energy(N);

    vector<vector<double>> s0 = determine_NonlinearConjugateGradient_s0(vertex); 
    vector<vector<double>> s1(vertex.size(), vector<double>(3,0.0)); 
    
    vector<int> timesOffNCG(N,0);
    bool IsCriteriaSatisfied = false; 
    bool updateReference = true;
    int iteration = 0;
    double TrialStepSize = 0.0;
    while ( IsCriteriaSatisfied == false && iteration < N-1){
        // The step size a0 needs a trial value, which is determined by rule-of-thumb. 
        if ( iteration == 0 || iteration % 50 == 0 ){
            vector<double> Scale_s0 = ForceScale(s0);
            double MaxForceScale = *max_element(Scale_s0.begin(), Scale_s0.end());
            TrialStepSize = l / MaxForceScale; // a0 = 1;
        }else{
            TrialStepSize = TrialStepSize * 2e1;
        } 
        {
            double StepSize = LineSearch_for_StepSize_to_minimize_energy(TrialStepSize, s0, vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
            cout<<"step: "<<iteration<<", trial StepSize = "<<TrialStepSize<<", StepSize = "<<StepSize<<endl;
            if ( StepSize == -1 ){
               cout<<"step: "<<iteration<<". Note: no efficent step size a is found. Stop now! "<<endl;
               //printoutREF(vertexref);
               //printoutstuck(vertex0);
               break;
            }
            update_vertex_from_NonlinearConjugateGradient_s0(vertex, face, StepSize, s0, param); 
            // calculate the new force and energy
            Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions);
            update_PreviousCoord_for_vertex(vertex); /// update the previous position
            update_PreviousForce_for_vertex(vertex);
            update_PreviousEnergy_for_face(face);
            param.EnergyPrevious = param.Energy;
            //////////////////////////////////////////////////////////////////////////////
            // calculate the direction s
            double NCGfactor0 = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor0) 
            for (int i = 0; i < vertex.size(); i++){
                vector<double> Forcetmp = - vertex[i].ForcePrevious.ForceTotal; 
                NCGfactor0 = NCGfactor0 + dot(Forcetmp, Forcetmp);
            }
            double NCGfactor = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor) 
            for (int i = 0; i < vertex.size(); i++){
                vector<double> Forcetmp = - vertex[i].Force.ForceTotal; 
                NCGfactor = NCGfactor + dot(Forcetmp, Forcetmp);
            }
            double peta1 = NCGfactor / NCGfactor0; if (param.usingNCG == false) peta1 = 0.0;
            #pragma omp parallel for 
            for (int i = 0; i < vertex.size(); i++){
                s1[i] = vertex[i].Force.ForceTotal + peta1 * s0[i];  
            }
            //////////////////////////////////////////////////////////////////////////
            // update the direction for Nonlinear Conjugate gradient method
            s0 = s1;
            TrialStepSize = StepSize;
            
            // store the energy and nodal force
            Energy[iteration] = param.Energy;
            MeanForce[iteration] = calculate_mean_force(vertex); 
            AreaTotal[iteration] = param.S;
            // update the reference 
            {
               if ( abs(Energy[iteration].EnergyTotal - Energy[iteration-1].EnergyTotal) < 1e-3 || abs(MeanForce[iteration]-MeanForce[iteration-1]) < 1e-3 ){
                   cout<<"update the reference structure! "<<endl;
                   update_reference_from_Coord(vertex);
               }
            }

            ////////////////////////////////////////////////////////
            // to check whether NCG has been stucked for too many times consistently. if so, change the flag to not-using-NCG 
            if ( param.isNCGstucked == true ){
                timesOffNCG[iteration] = 1;
            }
            if ( iteration > 10 ){
                int NoffNCG = accumulate(timesOffNCG.begin()+iteration-9, timesOffNCG.begin()+iteration, 0);
                if ( NoffNCG > 8 ){
                    param.usingNCG = false;
                }else if ( NoffNCG < 3 ){
                    param.usingNCG = true;
                    param.isNCGstucked = false;
                }
            }
        }
        // output parameters     
        cout<<"step: "<< iteration <<". Energy= "<<Energy[iteration].EnergyTotal<<". meanF= "<<MeanForce[iteration]<<endl;
        //////////////////////////////////////////////////////////////////////////////////
        // check whether to stop. if the total energy is flat for 100 simulation steps, then stop
        if ( iteration > 500 && abs((Energy[iteration].EnergyTotal-Energy[iteration-500].EnergyTotal)/500) < CriterionForEnergy && MeanForce[iteration] < CriterionForForce ){
            cout<<"The energy is minimized. Stop now!"<<endl;
            IsCriteriaSatisfied = true;
            break;
        }
        //double Ere_vs_Etot = Energy[iteration].EnergyRegularization / Energy[iteration].EnergyTolt;

        //////////////////////////////////////////////////////////////////////////////////
        // output the vertex
        if ( iteration % 100 == 0 ){    
            int kk = iteration/100;
            char filename[20] = "vertex%d.csv";
            sprintf(filename,"vertex%d.csv",kk);
            ofstream outfile(filename);
            for (int j = 0; j < vertex.size(); j++) {
                outfile << vertex[j].Coord[0] << ',' << vertex[j].Coord[1] << ',' << vertex[j].Coord[2] << '\n';
            }
            outfile.close();
        }
        // output energy and meanforce
        ofstream outfile2("EnergyForce.csv"); 
        for (int j = 0; j <= iteration; j++){
            if ( j == 0 ){
                outfile2 <<"Energy-Curvature, -Area, -Regularization, -Total ((pN.nm)); Mean Force (pN)"<< '\n';
            }
            outfile2 << Energy[j].EnergyCurvature << ", " << Energy[j].EnergyArea << ", " << Energy[j].EnergyRegularization << ", " << Energy[j].EnergyTotal<< ", " << MeanForce[j] << '\n';
        }
        outfile2.close();

        iteration ++;
    }
    // output the final structure
    ofstream outfile33("vertexfinal.csv");
    for (int j = 0; j < vertex.size(); j++) {
        outfile33 << setprecision(16) << vertex[j].CoordPrevious[0] << ',' << vertex[j].CoordPrevious[1] << ',' << vertex[j].CoordPrevious[2] << '\n';
    }
    outfile33.close();
}
