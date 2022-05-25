#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "vectorMath.cpp"
#pragma omp declare reduction(+: std::vector<double> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#pragma omp declare reduction(+: std::vector<vector<double>> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#pragma omp declare reduction(+: energy : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#pragma omp declare reduction(+: std::vector<force>  : omp_out += omp_in ) initializer( omp_priv = omp_orig )

using namespace std;

struct Vertex{
    int            Index;
    bool           IsOutLayer;
    bool           IsBoundary = false;
    bool           IsGhost = false;
    vector<double> Coord {0.0, 0.0, 0.0};
    vector<double> CoordPrevious {0.0, 0.0, 0.0};
    vector<double> ReferenceCoord {0.0, 0.0, 0.0};
    vector<int>    AdjacentFace; // faces that have this vertex. there should be 5 or 6 faces.
    vector<int>    AdjacentVertex; // vertice that are nearby the vertex_i. There should be 6 or less.
    force          ForcePrevious;
    force          Force;
};

struct Face{
    int         Index;
    bool        IsOutLayer;
    bool        IsBoundary = false;
    bool        IsGhost = false;
    bool        IsInsertionPatch = false;
    vector<int> AdjacentVertex {0, 0, 0};
    vector<int> OneRingVertex; // there should be 12 or 11 vertices.
    vector<int> AdjacentFace; // faces that are adjacent to this face. There should be 12 or 11 faces.
    double      SpontCurvature = 0.0;
    double      MeanCurvature = 0.0;
    vector<double> normVector {0.0, 0.0, 0.0}; // norm vector of this face element
    double      ElementArea = 0.0; // local area of this face.
    double      ElementVolume = 0.0; // local volume of this face.
    energy      EnergyPrevious;
    energy      Energy;
};

struct Deformation{
    int DeformShape = 0;
    int DeformArea  = 0;
    int Undeform    = 0;
};

struct Param {
    double kc;
    double us;
    double uv;
    double k;
    double K;
    double S0; //target area
    double S; // area, total area 
    double C0; // spontaneous curvature of insertions
    double c0; // spontaneous curvature of membrane
    double meanL;
    double gama_shape;
    double gama_area;
    double sigma = 0.0;
    int    GaussQuadratureN;
    bool   isInsertionAreaConstraint = false;
    bool   isAdditiveScheme = false;
    bool   isGlobalConstraint;
    double s0;
    bool   isBoundaryFixed;
    bool   isBoundaryPeriodic;
    bool   isBoundaryFree;
    double sideX;
    double sideY;
    double l;     // for flat subdivision 
    bool   usingNCG = true;
    bool   isNCGstucked = false;
    bool   usingRpi = true;
    vector<vector<int>> InsertionPatch;
    vector<int> AttachNode;
    vector<bool> isAttachNode;
    //mat    pointsposition;
    double bindCoefficient;
    double bindLength;
    double ClathRadius;
    double distToClath;
    vector<double> ClathCageCenter;

    energy Energy;
    energy EnergyPrevious;
    Deformation DeformationList;
}; 

vector<Vertex> setVertex_Loop_scheme(double sidex, double sidey, double l){ // vertex position
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    double lx = n * dx;
    double ly = m * dy;

    int nodenum = (n+1)*(m+1);
    vector<Vertex> vertex(nodenum);

    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        bool isEvenJ = false;
        if ( pow(-1.0,j) > 0.0 ){
            isEvenJ = true;
        }
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            double x = i*dx;
            if ( isEvenJ == true ){
                x = x + a/2.0;
            }
            double y = j*dy;
            vertex[index].Index = index;
            vertex[index].Coord[0] = x - lx/2.0; 
            vertex[index].Coord[1] = y - ly/2.0; 
            vertex[index].Coord[2] = 0.0;
        }
    }

    return vertex;
}
vector<Face> setFace_Loop_scheme(double sidex, double sidey, double l){ // face and its surrounding vertex
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int facenum = m*n*2;
    vector<Face> face(facenum);
    #pragma omp parallel for
    for (int j = 0; j < m; j++){
        bool isEvenJ = false;
        if ( pow(-1.0,j) > 0.0 ){
            isEvenJ = true;
        }
        for (int i = 0; i < n; i++){
            int index = 2*n*j + i*2;
            int node1, node2, node3, node4;
            if ( isEvenJ == false ){
                node1 = (n+1)*j + i;
                node2 = (n+1)*(j+1) + i;
                node3 = (n+1)*j + (i+1);
                node4 = (n+1)*(j+1) + (i+1);
            }else{
                node1 = (n+1)*(j+1) + i;
                node2 = (n+1)*(j+1) + (i+1);
                node3 = (n+1)*j + i;
                node4 = (n+1)*j + (i+1);
            }
            face[index].Index = index;
            face[index].AdjacentVertex[0] = node1;
            face[index].AdjacentVertex[1] = node2;
            face[index].AdjacentVertex[2] = node3;
            index = index + 1;
            face[index].Index = index;
            face[index].AdjacentVertex[0] = node2;
            face[index].AdjacentVertex[1] = node4;
            face[index].AdjacentVertex[2] = node3;
        }
    }
    return face;
}

void determine_Boundary_vertex_face(double sidex, double sidey, double l, vector<Vertex>& vertex, vector<Face>& face){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }

    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            if ( j == 0 || j == m || i == 0 || i == n ){
                vertex[index].IsBoundary = true; // if element is 1, then this vertex is on boundary
            }
        }
    }

    #pragma omp parallel for 
    for (int i = 0; i < face.size(); i++){
        int node1 = face[i].AdjacentVertex[0];
        int node2 = face[i].AdjacentVertex[1];
        int node3 = face[i].AdjacentVertex[2];
        if ( vertex[node1].IsBoundary == true || vertex[node2].IsBoundary == true || vertex[node3].IsBoundary == true ){
            face[i].IsBoundary = true; // this face is on boundary
        }
    }
}

void determine_Ghost_vertex_face(double sidex, double sidey, double l, bool isBoundaryFixed, bool isBoundaryPeriodic, bool isBoundaryFree, vector<Vertex>& vertex, vector<Face>& face){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    ///////////////////////////////////////// ghost vertex
    int vertexnum = (n+1)*(m+1);
    vector<int> TopBottom; 
    vector<int> LeftRight; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, 1, 2, m-2, m-1, m});
        LeftRight.insert(LeftRight.end(),{0, 1, 2, n-2, n-1, n});
    }else if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, m});
        LeftRight.insert(LeftRight.end(),{0, n});
    }
    // top and bottom ghost vertex
    for (int k = 0; k < TopBottom.size(); k++){
        int j = TopBottom[k];
        #pragma omp parallel for
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            vertex[index].IsGhost = true; 
        }
    }
    // left and right ghost vertex
    for (int k = 0; k < LeftRight.size(); k++){
        int i = LeftRight[k];
        #pragma omp parallel for
        for (int j = 0; j < m+1; j++){
            int index = (n+1)*j + i;
            vertex[index].IsGhost = true; 
        }
    }
    //////////////////////////////////////// ghost face
    int facenum = m*n*2;
    TopBottom.clear();; 
    LeftRight.clear();; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, 1, 2, m-3, m-2, m-1});
        LeftRight.insert(LeftRight.end(),{0, 1, 2, n-3, n-2, n-1});
    }else if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom.insert(TopBottom.end(),{0, m-1});
        LeftRight.insert(LeftRight.end(),{0, n-1});
    }
    for (int k = 0; k < TopBottom.size(); k++){
        int j = TopBottom[k];
        #pragma omp parallel for
        for (int i = 0; i < n; i++){
            int index = 2*n*j + i*2;
            face[index].IsGhost = true; 
            face[index+1].IsGhost = true;
        }
    }
    for (int k = 0; k < LeftRight.size(); k++){
        int i = LeftRight[k];
        #pragma omp parallel for
        for (int j = 0; j < m; j++){
            int index = 2*n*j + i*2;
            face[index].IsGhost = true; 
            face[index+1].IsGhost = true;
        }
    }
}

// find the faces around the vertex, probaly one, two, three or six faces that has vertex_i
void determine_AdjacentFace_for_vertex(vector<Vertex>& vertex, vector<Face>& face){
    #pragma omp parallel for 
     for (int i = 0; i < vertex.size(); i++){
        for (int j = 0; j < face.size(); j++){
            for (int k = 0; k < 3; k++){ // AdjacentVertex.size() = 3
                if ( i == face[j].AdjacentVertex[k] ){
                    vertex[i].AdjacentFace.push_back(j);
                }
            }
        }  
     }
}

// find out the nearby vertices around the vertex_i. there could be 6 or less.
void determine_AdjacentVertex_for_vertex(vector<Vertex>& vertex, vector<Face>& face){  
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        for (int j = 0; j < vertex[i].AdjacentFace.size(); j++){
            int FaceIndex = vertex[i].AdjacentFace[j]; 
            for (int k = 0; k < face[FaceIndex].AdjacentVertex.size(); k++){
                int VertexIndex = face[FaceIndex].AdjacentVertex[k];
                if ( VertexIndex != i ){
                    bool IsListed = false;
                    for (int m = 0; m < vertex[i].AdjacentVertex.size(); m++){
                        if ( VertexIndex == vertex[i].AdjacentVertex[m] ){
                            IsListed = true;
                        }
                    }
                    if ( IsListed == false ){
                        vertex[i].AdjacentVertex.push_back(VertexIndex);
                    }
                }
            }
        }
    }
}

int Find_NodeIndex(int node1, int node2, int node3, vector<Vertex>& vertex){
    int node = -1;
    for (int i = 0; i < vertex[node1].AdjacentVertex.size(); i++){
        int nodetmp1 = vertex[node1].AdjacentVertex[i];
        for (int j = 0; j < vertex[node2].AdjacentVertex.size(); j++){
            int nodetmp2 = vertex[node2].AdjacentVertex[j];
            if ( nodetmp1 == nodetmp2 && nodetmp1 != node3 ){
                node = nodetmp1;
            }
        }
    }
    if ( node == -1 ){
        cout<<"Wrong! No efficent NodeIndex is found in Find_NodeIndex. Node1 = "<<node1<<", Node2 = "<<node2<<", Node3 = "<<node3<<endl;
        exit(0);
    }
    return node;
}

// To find out the one-ring vertices aound face_i. It should be 12 for the flat surface because we set it up only with regular patch.
// The boundary faces do not have complete one-ring, neither it will be called in the code, so no need to store their one-ring-vertex
void determine_OneRingVertex_for_face(vector<Vertex>& vertex, vector<Face>& face){ 
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ){ 
            continue;
        }
        int d4 = face[i].AdjacentVertex[0]; 
        int d7 = face[i].AdjacentVertex[1]; 
        int d8 = face[i].AdjacentVertex[2]; 
        int d3  = Find_NodeIndex(d4, d7, d8, vertex);
        int d11 = Find_NodeIndex(d7, d8, d4, vertex);
        int d5  = Find_NodeIndex(d4, d8, d7, vertex);
        int d1  = Find_NodeIndex(d3, d4, d7, vertex);
        int d2  = Find_NodeIndex(d4, d5, d8, vertex);
        int d6  = Find_NodeIndex(d3, d7, d4, vertex);
        int d9  = Find_NodeIndex(d8, d5, d4, vertex);
        int d10 = Find_NodeIndex(d7, d11, d8, vertex);
        int d12 = Find_NodeIndex(d8, d11, d7, vertex);
        vector<int> tmp {d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12};
        face[i].OneRingVertex = tmp;
        //face[i].OneRingVertex.insert(face[i].OneRingVertex.end(),{d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12});
    }
}

void determine_IsInsertionPatch_for_face(vector<Face>& face, vector<vector<int>>& InsertionPatch){
    #pragma omp parallel for
    for (int i = 0; i < InsertionPatch.size(); i++){
        for (int j = 0; j < InsertionPatch[0].size(); j++){
            int FaceIndex = InsertionPatch[i][j];
            face[FaceIndex].IsInsertionPatch = true;
        }
    }
}

void determine_SpontaneousCurvature_for_face(double perturbed_c0, double unperturbed_c0, vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) 
            continue;

        if ( face[i].IsInsertionPatch == true ){
            face[i].SpontCurvature = perturbed_c0;
        }else{
            face[i].SpontCurvature = unperturbed_c0;
        }
    }
}

vector<vector<double>> setVMU(int n){ // To setup the Gaussian points for integral calculation. 'n' here is the Gausssian order.
    vector<vector<double>> vmu;
    if (n==1){
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0} };
        vmu = vmutmp;
    }else if(n==2){
        vector<vector<double> > vmutmp { {1.0/6.0, 1.0/6.0, 4.0/6.0},
                                         {1.0/6.0, 4.0/6.0, 1.0/6.0},
                                         {4.0/6.0, 1.0/6.0, 1.0/6.0} };
        vmu = vmutmp;
    }else if(n==3){ // third order Guasssian may not converge. Weird!
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0},
                                         {1.0/5.0, 1.0/5.0, 3.0/5.0},
                                         {1.0/5.0, 3.0/5.0, 1.0/5.0},
                                         {3.0/5.0, 1.0/5.0, 1.0/5.0} };
        /*
        vector<vector<double> > vmutmp { {1.0/3.0, 1.0/3.0, 1.0/3.0},
                                         {2.0/15.0, 11.0/15.0, 2.0/15.0},
                                         {2.0/15.0, 2.0/15.0, 11.0/15.0},
                                         {11.0/15.0, 2.0/15.0, 2.0/15.0} };
        */
        vmu = vmutmp;
    }else if(n==4){
        vector<vector<double> > vmutmp { {0.44594849091597, 0.44594849091597, 0.10810301816807},
                                         {0.44594849091597, 0.10810301816807, 0.44594849091597},
                                         {0.10810301816807, 0.44594849091597, 0.44594849091597},
                                         {0.09157621350977, 0.09157621350977, 0.81684757298046},
                                         {0.09157621350977, 0.81684757298046, 0.09157621350977},
                                         {0.81684757298046, 0.09157621350977, 0.09157621350977} };
        vmu = vmutmp;
    }else if(n==5){
        vector<vector<double> > vmutmp { {0.33333333333333, 0.33333333333333, 0.33333333333333},
                                         {0.47014206410511, 0.47014206410511, 0.05971587178977},
                                         {0.47014206410511, 0.05971587178977, 0.47014206410511},
                                         {0.05971587178977, 0.47014206410511, 0.47014206410511},
                                         {0.10128650732346, 0.10128650732346, 0.79742698535309},
                                         {0.10128650732346, 0.79742698535309, 0.10128650732346},
                                         {0.79742698535309, 0.10128650732346, 0.10128650732346} };
        vmu = vmutmp;
    }else if(n==6){
        vector<vector<double> > vmutmp { {0.24928674517091, 0.24928674517091, 0.50142650965818},
                                         {0.24928674517091, 0.50142650965818, 0.24928674517091},
                                         {0.50142650965818, 0.24928674517091, 0.24928674517091},
                                         {0.06308901449150, 0.06308901449150, 0.87382197101700},
                                         {0.06308901449150, 0.87382197101700, 0.06308901449150},
                                         {0.87382197101700, 0.06308901449150, 0.06308901449150},
                                         {0.31035245103378, 0.63650249912140, 0.05314504984482},
                                         {0.63650249912140, 0.05314504984482, 0.31035245103378},
                                         {0.05314504984482, 0.31035245103378, 0.63650249912140},
                                         {0.63650249912140, 0.31035245103378, 0.05314504984482},
                                         {0.31035245103378, 0.05314504984482, 0.63650249912140},
                                         {0.05314504984482, 0.63650249912140, 0.31035245103378} };
        vmu = vmutmp;
    }

    return vmu;
}

vector<double> setVMUcoefficient(int n){ // coefficeints for Gaussian points. 'n' here is the Gausssian order.
    vector<double> vmucoeff;
    if (n==1){
        vector<double> coefftmp { 1.0 };
        vmucoeff = coefftmp;
    }else if(n==2){
        vector<double> coefftmp { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
        vmucoeff = coefftmp;
    }else if(n==3){
        vector<double> coefftmp { -0.56250000000000, 0.52083333333333, 0.52083333333333, 0.52083333333333 };
        vmucoeff = coefftmp;
    }else if(n==4){
        vector<double> coefftmp { 0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532 };
        vmucoeff = coefftmp;
    }else if(n==5){
        vector<double> coefftmp { 0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483 };
        vmucoeff = coefftmp;
    }else if(n==6){
        vector<double> coefftmp { 0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021, 0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837 };
        vmucoeff = coefftmp;
    }

    return vmucoeff;
}

// sructure of shape funcitons.
struct Shapefunctions{
    vector<vector<double>> sf;
};
// shape functions on barycentric coordinates.
vector<vector<double>> determine_ShapeFunctions(vector<double>& vwu){
    double v = vwu[0]; double w = vwu[1]; double u = vwu[2];
    vector<vector<double>> sf(12,vector<double>(7)); // note, I transpose this matrix when output
    sf[0][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*v);
    sf[0][1] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf[0][2] = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf[0][3] = u*v; 
    sf[0][4] = pow(u,2.0) + u*v; 
    sf[0][5] = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf[0][6] = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf[1][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w); 
    sf[1][1] = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*w); 
    sf[1][2] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w);
    sf[1][3] = pow(u,2.0) + u*w; 
    sf[1][4] = u*w;
    sf[1][5] = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf[1][6] = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf[2][0] = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w + 6.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 12.0*pow(u,2.0)*pow(v,2.0) + 6.0*u*pow(v,2.0)*w + 6.0*u*pow(v,3.0) + 2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[2][1] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*v - 6.0*u*pow(v,2.0) - 2.0*pow(v,3.0));
    sf[2][2] = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w - 12.0*pow(u,2)*v - 12.0*u*v*w - 18.0*u*pow(v,2.0) - 6.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf[2][3] = -2.0*u*v;
    sf[2][4] = u*w + v*w + u*v + pow(v,2.0);
    sf[2][5] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf[2][6] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf[3][0] = 1.0/12.0*(6.0*pow(u,4.0) + 24.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 8.0*u*pow(w,3.0) + pow(w,4.0) + 24.0*pow(u,3.0)*v + 60.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 8.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[3][1] = 1.0/12.0*(-12.0*pow(u,2.0)*w - 12.0*u*pow(w,2.0) - 2.0*pow(w,3.0) - 24.0*pow(u,2.0)*v - 48.0*u*v*w - 12.0*v*pow(w,2.0) -24.0*u*pow(v,2.0) - 18.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf[3][2] = 1.0/12.0*(-24.0*pow(u,2.0)*w - 24.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 12.0*pow(u,2.0)*v - 48.0*u*v*w - 18.0*v*pow(w,2.0) - 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
    sf[3][3] = -2.0*u*w - 2.0*pow(u,2.0) + v*w + pow(v,2.0);
    sf[3][4] = -2.0*pow(u,2.0) + pow(w,2.0) - 2.0*u*v + v*w;
    sf[3][5] = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf[3][6] = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf[4][0] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0) + pow(w,4.0) + 2.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 6.0*u*v*pow(w,2.0) + 2.0*v*pow(w,3.0));
    sf[4][1] = 1.0/12.0*(-2.0*pow(u,3.0) - 12.0*pow(u,2.0)*w - 18.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 6.0*pow(u,2.0)*v - 12.0*u*v*w - 6.0*v*pow(w,2.0));
    sf[4][2] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0));
    sf[4][3] = u*w + pow(w,2.0) + u*v + v*w;
    sf[4][4] = -2.0*u*w;
    sf[4][5] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf[4][6] = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf[5][0] = 1.0/12.0*(2.0*u*pow(v,3.0) + pow(v,4.0)); 
    sf[5][1] = 1.0/12.0*(6.0*u*pow(v,2.0) + 2.0*pow(v,3.0)); 
    sf[5][2] = -1.0/6.0*pow(v,3.0);
    sf[5][3] = u*v; 
    sf[5][4] = 0.0;
    sf[5][5] = -1.0/2.0*pow(v,2.0); 
    sf[5][6] = -1.0/2.0*pow(v,2.0);
    sf[6][0] = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0)+ pow(w,4.0) + 8.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 8.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 60.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 24.0*u*pow(v,3.0) + 24.0*pow(v,3.0)*w + 6.0*pow(v,4.0));
    sf[6][1] = 1.0/12.0*(4.0*pow(u,3.0) + 18.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 2.0*pow(w,3.0) + 24.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 24.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w);
    sf[6][2] = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0) + 12.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w);
    sf[6][3] = pow(u,2.0) + u*w - 2.0*v*w - 2.0*pow(v,2.0);
    sf[6][4] = -2.0*u*w - 2.0*u*v - 2.0*v*w - 2.0*pow(v,2.0);
    sf[6][5] = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
    sf[6][6] = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
    sf[7][0] = 1.0/12.0*(pow(u,4.0) + 8.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 24.0*u*pow(w,3.0) + 6.0*pow(w,4.0) + 6.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 60.0*u*v*pow(w,2.0) + 24.0*v*pow(w,3.0) + 12.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 6.0*u*pow(v,3.0) + 8.0*pow(v,3.0)*w + pow(v,4.0));
    sf[7][1] = 1.0/12.0*(2.0*pow(u,3.0) + 12.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 6.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) - 6.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
    sf[7][2] = 1.0/12.0*(4.0*pow(u,3.0) + 24.0*pow(u,2.0)*w + 24.0*u*pow(w,2.0) + 18.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w + 2.0*pow(v,3.0));
    sf[7][3] = -2.0*u*w - 2.0*pow(w,2.0) - 2.0*u*v - 2.0*v*w;
    sf[7][4] = pow(u,2.0) - 2.0*pow(w,2.0) + u*v - 2.0*v*w;
    sf[7][5] = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
    sf[7][6] = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
    sf[8][0] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0)); 
    sf[8][1] = -1.0/6.0*pow(w,3.0); 
    sf[8][2] = 1.0/12.0*(6.0*u*pow(w,2.0) + 2.0*pow(w,3.0));
    sf[8][3] = 0.0; 
    sf[8][4] = u*w;
    sf[8][5] = -1.0/2.0*pow(w,2.0); 
    sf[8][6] = -1.0/2.0*pow(w,2.0);
    sf[9][0] = 1.0/12.0*(2.0*pow(v,3.0)*w + pow(v,4.0));
    sf[9][1] = 1.0/12.0*(6.0*pow(v,2.0)*w + 4.0*pow(v,3.0)); 
    sf[9][2] = 1.0/6.0*pow(v,3.0);
    sf[9][3] = v*w + pow(v,2.0); 
    sf[9][4] = 0.0;
    sf[9][5] = 1.0/2.0*pow(v,2.0); 
    sf[9][6] = 1.0/2.0*pow(v,2.0);
    sf[10][0] = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0) + 6.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 6.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 2.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf[10][1] = 1.0/12.0*(4.0*pow(w,3.0) + 18.0*v*pow(w,2.0) + 6.0*u*pow(w,2.0) + 12.0*pow(v,2.0)*w + 12.0*u*v*w + 2.0*pow(v,3.0) + 6.0*u*pow(v,2.0));
    sf[10][2] = 1.0/12.0*(2.0*pow(w,3.0) + 6.0*u*pow(w,2.0) + 12.0*v*pow(w,2.0) + 12.0*u*v*w + 18.0*pow(v,2.0)*w + 6.0*u*pow(v,2.0) + 4.0*pow(v,3.0));
    sf[10][3] = pow(w,2.0) + v*w + u*w + u*v;
    sf[10][4] = u*w + v*w + u*v + pow(v,2.0);
    sf[10][5] = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
    sf[10][6] = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
    sf[11][0] = 1.0/12.0*(pow(w,4.0) + 2.0*v*pow(w,3.0)); 
    sf[11][1] = 1.0/6.0*pow(w,3.0); 
    sf[11][2] = 1.0/12.0*(4.0*pow(w,3.0) + 6.0*v*pow(w,2.0));
    sf[11][3] = 0.0; 
    sf[11][4] = pow(w,2.0) + v*w;
    sf[11][5] = 1.0/2.0*pow(w,2.0); 
    sf[11][6] = 1.0/2.0*pow(w,2.0);

    return transpose(sf);
    // 12 shape functions and their differential equations; 
    // shape_functions(1,:), shape functions;  
    // shape_functions(2,:), differential to v; 
    // shape_functions(3,:), differential to w; 
    // shape_functions(4,:), double differential to v; 
    // shape_functions(5,:), double differential to w;
    // shape_functions(6,:), differential to v and w; 
    // shape_functions(7,:), differential to w and v;
}
/*
// for irregular patch, more subdivision is needed. 
// for different sub-element, different new nodes are selected, select-matrix (SM)
// vertex here is the original 11 vertice, so vertex is 11*3 matrix
// M(17,11); mat M1(12,17); mat M2(12,17); mat M3(12,17); mat M4(11,17);
void determine_SubdivisionMatrix(Matrix& M, Matrix& SM1, Matrix& SM2, Matrix& SM3, Matrix& SM4){
    int N = 6; double w = 3.0/8.0/N; // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
    int N1 = 5; double w1 = 3.0/8.0/N1; // w1=1/N1*(5/8-(3/8+1/4*cos(2*pi/N1))^2);
    double a = 3.0/8.0; double b = 1.0/8.0;
    vector<vector<double> > Mtmp { {a, b, a, b, 0, 0, 0, 0, 0, 0, 0},
                                   {b, a, a, 0, 0, b, 0, 0, 0, 0, 0},
                                   {w1, w1, 1.0-N1*w1, w1, 0, w1, w1, 0, 0, 0, 0},
                                   {b, 0, a, a, 0, 0, b, 0, 0, 0, 0},
                                   {0, a, b, 0, b, a, 0, 0, 0, 0, 0},
                                   {0, b, a, 0, 0, a, b, 0, 0, 0, 0},
                                   {0, 0, a, b, 0, b, a, 0, 0, 0, 0},
                                   {0, 0, b, a, 0, 0, a, b, 0, 0, 0},
                                   {0, b, 0, 0, a, a, 0, 0, b, 0, 0},
                                   {0, w, w, 0, w, 1.0-N*w,  w, 0, w, w, 0},
                                   {0, 0, b, 0, 0, a, a, 0, 0, b, 0},
                                   {0, 0, w, w, 0, w, 1.0-N*w, w, 0, w, w},
                                   {0, 0, 0, b, 0, 0, a, a, 0, 0, b},
                                   {0, 0, 0, 0, b, a, 0, 0, a, b, 0},
                                   {0, 0, 0, 0, 0, a, b, 0, b, a, 0},
                                   {0, 0, 0, 0, 0, b, a, 0, 0, a, b},
                                   {0, 0, 0, 0, 0, 0, a, b, 0, b, a} };
    M = Mtmp;
    // SM1=zeros(12,17);
    vector<int> element1 { 2, 3, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16 };
    for (int i = 0; i < 12; i++){
        SM1(i,element1[i]) = 1.0;
    }
    // SM2=zeros(12,17);
    vector<int> element2 { 4, 1, 9, 5, 2, 14, 10, 6, 3, 15, 11, 7 };
    for (int i = 0; i < 12; i++){
        SM2(i,element2[i]) = 1.0;
    }
    // SM3=zeros(12,17);
    vector<int> element3 { 1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15 };
    for (int i = 0; i < 12; i++){ 
        SM3(i,element3[i]) = 1.0;
    }
    // SM4=zeros(11,17);
    vector<int> element4 { 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11 };
    for (int i = 0; i < 11; i++){ 
        SM4(i,element4[i]) = 1.0;
    }
}
*/
void calculate_element_area_volume(vector<Vertex>& vertex, vector<Face>& face, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions){ 
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        
        int numberOneRingVertex = face[i].OneRingVertex.size(); // 12
        vector<vector<double>> dots(numberOneRingVertex,vector<double>(3)); //dots(12,3); 12 nodes
        for (int j = 0; j < numberOneRingVertex; j++){
            int NodeIndex = face[i].OneRingVertex[j];
            dots[j] = vertex[NodeIndex].Coord;
        }
        // Gaussian quadrature, 3 points
        double area = 0.0;
        double volume = 0.0;
        for (int j = 0; j < GaussQuadratureCoeff.size(); j++){
            vector<vector<double>> sf = ShapeFunctions[j].sf;
            vector<double> x = sf[0] * dots;
            vector<double> a_1 = sf[1] * dots;
            vector<double> a_2 = sf[2] * dots;
            vector<double> a_3 = cross(a_1,a_2);
            double  sqa = norm(a_3); 
            double s = sqa; 
            vector<double> d = a_3 / sqa;
            area = area + 1.0/2.0 * GaussQuadratureCoeff[j] * s; 

            double v = s * dot(x,d);
            volume = volume + 1.0/2.0 * GaussQuadratureCoeff[j] * v; 
        }

        face[i].ElementArea = area;
        face[i].ElementVolume = volume;
    }
}

double sum_Membrane_Area(vector<Face>& face){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        sum += face[i].ElementArea;
    }
    return sum;
}
double sum_Membrane_Volume(vector<Face>& face){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        sum += face[i].ElementVolume;
    }
    return sum;
}
void clear_forceONvertex_and_energyONface(vector<Vertex>& vertex, vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        force Forcetmp;
        vertex[i].Force = Forcetmp;
    }
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        energy Energytmp;
        face[i].Energy = Energytmp;
    }
}
/*
void calculate_total_energy_save_in_paramEnergy(vector<Face>& face, Param& param){
    energy tmpEnergy;
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        tmpEnergy.EnergyCurvature      = tmpEnergy.EnergyCurvature + face[i].Energy.EnergyCurvature;
        tmpEnergy.EnergyArea           = tmpEnergy.EnergyArea + face[i].Energy.EnergyArea;
        tmpEnergy.EnergyVolume         = tmpEnergy.EnergyVolume + face[i].Energy.EnergyVolume;
        tmpEnergy.EnergyThickness      = tmpEnergy.EnergyThickness + face[i].Energy.EnergyThickness;
        tmpEnergy.EnergyTilt           = tmpEnergy.EnergyTilt + face[i].Energy.EnergyTilt;
        tmpEnergy.EnergyRegularization = tmpEnergy.EnergyRegularization + face[i].Energy.EnergyRegularization;
    }
    double EnergyArea_GlobalConstraint = 0.0;
    double us = param.us;
    double S0 = param.S0;
    double S  = param.S;
    EnergyArea_GlobalConstraint = 0.5*us/S0*pow(S-S0,2.0); 
    if (param.isGlobalConstraint == true)
        tmpEnergy.EnergyArea = EnergyArea_GlobalConstraint;

    param.Energy.EnergyCurvature      = tmpEnergy.EnergyCurvature;
    param.Energy.EnergyArea           = tmpEnergy.EnergyArea; 
    param.Energy.EnergyVolume         = tmpEnergy.EnergyVolume;
    param.Energy.EnergyThickness      = tmpEnergy.EnergyThickness;
    param.Energy.EnergyTilt           = tmpEnergy.EnergyTilt;
    param.Energy.EnergyRegularization = tmpEnergy.EnergyRegularization;
    param.Energy.EnergyTotal          = tmpEnergy.EnergyCurvature + tmpEnergy.EnergyArea + tmpEnergy.EnergyVolume + tmpEnergy.EnergyThickness + tmpEnergy.EnergyTilt + tmpEnergy.EnergyRegularization;
}
*/
void update_PreviousCoord_for_vertex(vector<Vertex>& vertex){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].CoordPrevious = vertex[i].Coord;
    }
}

void update_reference_from_Coord(vector<Vertex>& vertex){
    #pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].ReferenceCoord = vertex[i].CoordPrevious;
    }
}

void update_PreviousForce_for_vertex(vector<Vertex>& vertex){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].ForcePrevious = vertex[i].Force;
    }
}

void update_PreviousEnergy_for_face(vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        face[i].EnergyPrevious = face[i].Energy;
    }
}
vector<vector<double>> determine_NonlinearConjugateGradient_s0(vector<Vertex>& vertex){
    vector<vector<double>> s0(vertex.size(), vector<double>(3,0.0));
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        s0[i] = vertex[i].Force.ForceTotal;
    }
    return s0;
}

vector<double> ForceScale(vector<vector<double>>& force){
    int num = force.size();
    vector<double> out(num);
    #pragma omp parallel for
    for (int i = 0; i < num; i++){
        out[i] = norm( force[i]);
    }
    return out;
}

double calculate_mean_force(vector<Vertex>& vertex){
    vector<double> forcescale(vertex.size(),0.0);
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        forcescale[i] = norm( vertex[i].Force.ForceTotal );
    }
    double sum = 0.0;
    #pragma omp parallel for reduction (+:sum)
    for (int i = 0; i < vertex.size(); i++){
        sum += forcescale[i];
    }
    return sum / vertex.size();
}


void manage_force_for_boundary_ghost_vertex(vector<Vertex>& vertex, vector<Face>& face, Param& param){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;
    if ( param.isBoundaryFixed == true ){
        force zeroForce;
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( face[i].IsBoundary == false )
                continue;
            int node1 = face[i].AdjacentVertex[0];
            int node2 = face[i].AdjacentVertex[1];
            int node3 = face[i].AdjacentVertex[2];
            vertex[node1].Force = zeroForce;
            vertex[node2].Force = zeroForce;
            vertex[node3].Force = zeroForce;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        force zeroForce;
        #pragma omp parallel for 
        for (int i = 0; i < vertexnum; i++){
            if ( vertex[i].IsGhost == true ){
                vertex[i].Force = zeroForce;
            }
        }
    }else if ( param.isBoundaryFree == true ){   
        force zeroForce;
        #pragma omp parallel for 
        for ( int i = 0; i < n+1; i++ ){
            int j = 0;
            int index = (n+1)*j + i;
            vertex[index].Force = zeroForce;
            index = (n+1)*m + i;
            vertex[index].Force = zeroForce;
        }
        #pragma omp parallel for
        for ( int j = 0; j < m+1; j++ ){
            int i = 0;
            int index = (n+1)*j + i;
            vertex[index].Force = zeroForce;
            index = (n+1)*j + n;
            vertex[index].Force = zeroForce;
        }
    }
}


void update_vertex_from_NonlinearConjugateGradient_s0(vector<Vertex>& vertex, vector<Face>& face, double a, vector<vector<double>>& Force, Param& param){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double aa = dx;
    double dy = sqrt(3.0)/2.0 * aa; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;

    // update the vertex position Coord.
    #pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].Coord = vertex[i].CoordPrevious + a * Force[i];
    }
    // deal with the boundar/ghost vertex
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( face[i].IsBoundary != true )
                continue;
            int node1 = face[i].AdjacentVertex[0];
            int node2 = face[i].AdjacentVertex[1];
            int node3 = face[i].AdjacentVertex[2];
            vertex[node1].Coord = vertex[node1].CoordPrevious;
            vertex[node2].Coord = vertex[node2].CoordPrevious;
            vertex[node3].Coord = vertex[node3].CoordPrevious;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index0 = n*0 + i;
            int index1 = n*1 + i;
            int index2 = n*2 + i;
            int index00 = n*(m-6) + i;
            int index11 = n*(m-5) + i;
            int index22 = n*(m-4) + i;
            vertex[index0].Coord = vertex[index0].CoordPrevious + ( vertex[index00].Coord-vertex[index00].CoordPrevious );
            vertex[index1].Coord = vertex[index1].CoordPrevious + ( vertex[index11].Coord-vertex[index11].CoordPrevious );
            vertex[index2].Coord = vertex[index2].CoordPrevious + ( vertex[index22].Coord-vertex[index22].CoordPrevious );
        }
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index0 = n*(m-2) + i;
            int index1 = n*(m-1) + i;
            int index2 = n*(m-0) + i;
            int index00 = n*4 + i;
            int index11 = n*5 + i;
            int index22 = n*6 + i;
            vertex[index0].Coord = vertex[index0].CoordPrevious + ( vertex[index00].Coord-vertex[index00].CoordPrevious );
            vertex[index1].Coord = vertex[index1].CoordPrevious + ( vertex[index11].Coord-vertex[index11].CoordPrevious );
            vertex[index2].Coord = vertex[index2].CoordPrevious + ( vertex[index22].Coord-vertex[index22].CoordPrevious );
        }
    }else if ( param.isBoundaryFree == true ){
        int index1, index2, index3, index4;
        // left side
        for (int j = 2; j < m ; j++ ){
            if ( pow(-1.0,j) > 0.0 ){
                index1 = (n+1)*j + 0;
                index2 = (n+1)*j + 1;
                index3 = (n+1)*(j-1) + 1;
                index4 = (n+1)*(j-1) + 2;
            }else{
                index1 = (n+1)*j + 0;
                index2 = (n+1)*j + 1;
                index3 = (n+1)*(j-1) + 0;
                index4 = (n+1)*(j-1) + 1;
            }
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*1 + 0; 
        index2 = (n+1)*2 + 0;
        index3 = (n+1)*1 + 1;
        index4 = (n+1)*2 + 1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        // right side
        index1 = (n+1)*1 + n; 
        index2 = (n+1)*2 + n-1;
        index3 = (n+1)*1 + n-1;
        index4 = (n+1)*2 + n-2;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        for (int j = 2; j < m ; j++ ){
            if ( pow(-1.0,j) > 0.0 ){
                index1 = (n+1)*j + n;
                index2 = (n+1)*j + n-1;
                index3 = (n+1)*(j-1) + n;
                index4 = (n+1)*(j-1) + n-1;
            }else{
                index1 = (n+1)*j + n;
                index2 = (n+1)*j + n-1;
                index3 = (n+1)*(j-1) + n-1;
                index4 = (n+1)*(j-1) + n-2;
            }
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        // bottom 
        for (int i = 0; i < n; i++){
            index1 = (n+1)*0 + i;
            index2 = (n+1)*1 + i;
            index3 = (n+1)*1 + i+1;
            index4 = (n+1)*2 + i;
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*0 + n; 
        index2 = (n+1)*0 + n-1;
        index3 = (n+1)*1 + n;
        index4 = (n+1)*1 + n-1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        // top 
        for (int i = 1; i < n+1; i++){
            index1 = (n+1)*m + i;
            index2 = (n+1)*(m-1) + i;
            index3 = (n+1)*(m-1) + i+1;
            index4 = (n+1)*(m-2) + i;
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*m + n; 
        index2 = (n+1)*m + n-1;
        index3 = (n+1)*(m-1) + n;
        index4 = (n+1)*(m-1) + n-1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
    } 
}
/*
struct Particle{
    int Index;
    vector<double> Coord {0.0, 0.0, 0.0};
    int FaceIndex; // face index that this particle is located on
    double D; // diffusion constant; nm2/us
}
void single_particles_diffusing_on_surface(vector<Vertex>& vertex, vector<Face>& face, Param& param){
    int CopyNumber = 1000;
    vector<double> TitrationCenter {param.sideX, param.sideY, 0.0};
    int TimeSteps = 1E5; // us
    double dt =  0.1; // us
    double D = 1.0;  // nm2/us
    // initialize all the particles
    for (int i = 0; i < Copy)
    for (int i = 0; i < TimeSteps; i++){
        #pragma omp parallel for
        for(int j = 0; j < CopyNumber; j++){

        }
    }

}
*/