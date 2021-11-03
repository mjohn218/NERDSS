#include <math.h>
#include <cmath>
#include <iostream>
#include <armadillo>
#include <omp.h>
#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) initializer( omp_priv = omp_orig )

using namespace std;
using namespace arma;

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
    double s0;
    bool   isBoundaryFixed;
    bool   isBoundaryPeriodic;
    bool   isBoundaryFree;
    double sideX;
    double sideY;
    double l;     // for flat subdivision 
    bool   usingNCG;
    bool   usingRpi;
    Mat<int> insertionpatch;
    Row<int> AttachNode;
    vector<bool> isAttachNode;
    mat    pointsposition;
    double bindCoefficient;
    double bindLength;
    double ClathRadius;
    double distToClath;
    vector<double> ClathCageCenter;
}; 

Mat<double> setvertex_Loop_scheme(double sidex, double sidey, double l){ // vertex position
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    double lx = n * dx;
    double ly = m * dy;

    int nodenum = (n+1)*(m+1);
    mat vertex(nodenum,3); vertex.fill(0.0);
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
            vertex(index,0) = x - lx/2.0; 
            vertex(index,1) = y - ly/2.0;
        }
    }
    return vertex;
}
Mat<int> setface_Loop_scheme(double sidex, double sidey, double l){ // face and its surrounding vertex
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int facenum = m*n*2;
    Mat<int> face(facenum,3);
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
            face(index,0) = node1; face(index,1) = node2; face(index,2) = node3;
            face(index+1,0) = node2; face(index+1,1) = node4; face(index+1,2) = node3;
        }
    }
    return face;
}

Row<int> determine_BoundaryVertex(double sidex, double sidey, double l){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    Row<int> isBoundaryNode(vertexnum); isBoundaryNode.fill(-1);
    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            if ( j == 0 || j == m || i == 0 || i == n ){
                isBoundaryNode(index) = 1; // if element is 1, then this vertex is on boundary
            }
        }
    }
    return isBoundaryNode;
}

Row<int> determine_Boundaryface(Mat<int> face, Row<int> isBoundaryNode){
    int facenum = face.n_rows;
    Row<int> isBoundaryFace(facenum); isBoundaryFace.fill(-1);
    #pragma omp parallel for 
    for (int i = 0; i < facenum; i++){
        int node1 = face(i,0);
        int node2 = face(i,1);
        int node3 = face(i,2);
        bool isboundaryface = false;
        if ( isBoundaryNode(node1) == 1 || isBoundaryNode(node2) == 1 || isBoundaryNode(node3) == 1 ){
            isboundaryface = true;
            isBoundaryFace(i) = 1; // if element is 1, then this face is on boundary
        }
    }
    return isBoundaryFace;
}

Row<int> determine_GhostVertex(double sidex, double sidey, double l, bool isBoundaryFixed, bool isBoundaryPeriodic, bool isBoundaryFree){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    Row<int> isGhostNode(vertexnum); isGhostNode.fill(-1);

    Row<int> TopBottom; 
    Row<int> LeftRight; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom << 0 << 1 << 2 << m-2 << m-1 << m << endr;
        LeftRight << 0 << 1 << 2 << n-2 << n-1 << n << endr;
    }
    if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom << 0 << m << endr;
        LeftRight << 0 << n << endr;
    }
    int number = TopBottom.n_cols;
    for (int k = 0; k < number; k++){
        int j = TopBottom(k);
        #pragma omp parallel for
        for (int i = 0; i < n+1; i++){
            int index = (n+1)*j + i;
            isGhostNode(index) = 1; // if element is 1, then this vertex is ghost
        }
    }
    for (int k = 0; k < number; k++){
        int i = LeftRight(k);
        #pragma omp parallel for
        for (int j = 0; j < m+1; j++){
            int index = (n+1)*j + i;
            isGhostNode(index) = 1; // if element is 1, then this vertex is ghost
        }
    }
    return isGhostNode;
}

Row<int> determine_GhostFace(double sidex, double sidey, double l, bool isBoundaryFixed, bool isBoundaryPeriodic, bool isBoundaryFree){
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int facenum = m*n*2;
    Row<int> isGhostFace(facenum); isGhostFace.fill(-1);

    Row<int> TopBottom; 
    Row<int> LeftRight; 
    if (isBoundaryPeriodic == true && isBoundaryFree == false && isBoundaryFixed == false){
        TopBottom << 0 << 1 << 2 << m-3 << m-2 << m-1 << endr;
        LeftRight << 0 << 1 << 2 << n-3 << n-2 << n-1 << endr;
    }
    if (isBoundaryFree == true && isBoundaryPeriodic == false && isBoundaryFixed == false){
        TopBottom << 0 << m-1 << endr;
        LeftRight << 0 << n-1 << endr;
    }
    int number = TopBottom.n_cols;
    for (int k = 0; k < number; k++){
        int j = TopBottom(k);
        #pragma omp parallel for
        for (int i = 0; i < n; i++){
            int index = 2*n*j + i*2;
            isGhostFace(index) = 1; // if element is 1, then this face is ghost
            isGhostFace(index+1) = 1;
        }
    }
    for (int k = 0; k < number; k++){
        int i = LeftRight(k);
        #pragma omp parallel for
        for (int j = 0; j < m; j++){
            int index = 2*n*j + i*2;
            isGhostFace(index) = 1; // if element is 1, then this face is ghost
            isGhostFace(index+1) = 1;
        }
    }

    return isGhostFace;
}

// find the faces around one vertex, probaly one, two, three or six faces that has vertex_i
Mat<int> vertexi_face_with_it(mat vertex, Mat<int> face){
    int vertexnum = vertex.n_rows;
    Mat<int> vertexi_face(vertexnum,6); vertexi_face.fill(-1);
    #pragma omp parallel for 
     for (int i = 0; i < vertex.n_rows; i++){
        int facenumber = -1;
        for (int j = 0; j < face.n_rows; j++){
            for (int k = 0; k < face.n_cols; k++){
                if ( i == face(j,k) ){
                    facenumber = facenumber + 1;
                    vertexi_face(i,facenumber) = j;
                }
            }
        }  
     }
     return vertexi_face;
}

// find out the nearby vertices around the vertex_i. there could be 6 or less.
Mat<int> neighbor_vertices(Mat<int> vertexi_faces, Mat<int> face){  
    int vertexnum = vertexi_faces.n_rows;
    Mat<int> nearby_vertices(vertexnum,6); nearby_vertices.fill(-1);
    #pragma omp parallel for
    for (int i = 0; i < vertexnum; i++){
        int N = 6;
        Row<int> A(N); A.fill(-1);
        int shu = -1;
        for (int j = 0; j < N; j++){
            int Numface = vertexi_faces(i,j); 
            if (Numface == -1)
                continue;  
            
            for (int k = 0; k < 3; k++){
                int Numvertex = face(Numface,k);
                if ( Numvertex != i ){
                    bool islisted = false;
                    for (int m = 0; m < N; m++){
                        if ( Numvertex == A(m) ){
                            islisted = true;
                        }
                    }
                    if ( islisted == false ){
                         shu = shu + 1;
                         A(shu) = Numvertex;
                    }
                }
            }
        }
        for (int j = 0; j < N; j++){
            nearby_vertices(i,j) = A(j);
        }
    }
    return nearby_vertices;
}

int find_nodeindex(int node1, int node2, int node3, Mat<int> closest_nodes){
    int node;
    for (int i = 0; i < closest_nodes.n_cols; i++){
        int nodetmp1 = closest_nodes(node1,i);
        for (int j = 0; j < closest_nodes.n_cols; j++){
            int nodetmp2 = closest_nodes(node2,j);
            if ( nodetmp1 == nodetmp2 && nodetmp1 != -1 && nodetmp1 != node3 ){
                node = nodetmp1;
            }
        }
    }
    return node;
}

// find out the one-ring vertices aound face_i. It could be 12 or less. 
Mat<int> one_ring_vertices(Mat<int> face, Mat<int> closest_nodes, Row<int> isBoundaryFace){ 
    int facenum = face.n_rows;
    Mat<int>  ring_vertices(facenum,12);  ring_vertices.fill(-1);
    #pragma omp parallel for
    for (int i = 0; i < facenum; i++){
        if ( isBoundaryFace(i) == 1 ){
            continue;
        }
        int d4 = face(i,0); 
        int d7 = face(i,1);
        int d8 = face(i,2);
        int d3 = find_nodeindex(d4, d7, d8, closest_nodes);
        int d11 = find_nodeindex(d7, d8, d4, closest_nodes);
        int d5 = find_nodeindex(d4, d8, d7, closest_nodes);
        int d1 = find_nodeindex(d3, d4, d7, closest_nodes);
        int d2 = find_nodeindex(d4, d5, d8, closest_nodes);
        int d6 = find_nodeindex(d3, d7, d4, closest_nodes);
        int d9 = find_nodeindex(d8, d5, d4, closest_nodes);
        int d10 = find_nodeindex(d7, d11, d8, closest_nodes);
        int d12 = find_nodeindex(d8, d11, d7, closest_nodes);
        Row<int> v(12); v << d1 << d2 << d3 << d4 << d5 << d6 << d7 << d8 << d9 << d10 << d11 << d12;
        ring_vertices.row(i) = v;
    }
    return ring_vertices;
}

void shapefunctions(rowvec vwu, mat& sf){
    // 12 shape functions and their differential equations; shape_functions(:,1), shape functions;  
    // shape_functions(:,2), differential to v; shape_functions(:,3), differential to w; 
    // shape_functions(:,4), double differential to v; shape_functions(:,5), double differential to w;
    // shape_functions(:,6), differential to v and w; shape_functions(:,7), differential to w and v;
    double v = vwu(0); double w = vwu(1); double u = vwu(2);
    //sf=zeros(12,7);
    sf(0,0) = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*v);
    sf(0,1) = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf(0,2) = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*v); 
    sf(0,3) = u*v; 
    sf(0,4) = pow(u,2.0) + u*v; 
    sf(0,5) = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf(0,6) = 1.0/2.0*(pow(u,2.0) + 2.0*u*v); 
    sf(1,0) = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w); 
    sf(1,1) = 1.0/12.0*(-4.0*pow(u,3.0) - 6.0*pow(u,2.0)*w); 
    sf(1,2) = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w);
    sf(1,3) = pow(u,2.0) + u*w; 
    sf(1,4) = u*w;
    sf(1,5) = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf(1,6) = 1.0/2.0*(pow(u,2.0) + 2.0*u*w); 
    sf(2,0) = 1.0/12.0*(pow(u,4.0) + 2.0*pow(u,3.0)*w + 6.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 12.0*pow(u,2.0)*pow(v,2.0) + 6.0*u*pow(v,2.0)*w + 6.0*u*pow(v,3.0) + 2.0*pow(v,3.0)*w + pow(v,4.0));
    sf(2,1) = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*v - 6.0*u*pow(v,2.0) - 2.0*pow(v,3.0));
    sf(2,2) = 1.0/12.0*(-2.0*pow(u,3.0) - 6.0*pow(u,2.0)*w - 12.0*pow(u,2)*v - 12.0*u*v*w - 18.0*u*pow(v,2.0) - 6.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf(2,3) = -2.0*u*v;
    sf(2,4) = u*w + v*w + u*v + pow(v,2.0);
    sf(2,5) = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf(2,6) = 1.0/2.0*(-pow(u,2.0) - 2.0*u*v + pow(v,2.0));
    sf(3,0) = 1.0/12.0*(6.0*pow(u,4.0) + 24.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 8.0*u*pow(w,3.0) + pow(w,4.0) + 24.0*pow(u,3.0)*v + 60.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 8.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
    sf(3,1) = 1.0/12.0*(-12.0*pow(u,2.0)*w - 12.0*u*pow(w,2.0) - 2.0*pow(w,3.0) - 24.0*pow(u,2.0)*v - 48.0*u*v*w - 12.0*v*pow(w,2.0) -24.0*u*pow(v,2.0) - 18.0*pow(v,2.0)*w - 4.0*pow(v,3.0));
    sf(3,2) = 1.0/12.0*(-24.0*pow(u,2.0)*w - 24.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 12.0*pow(u,2.0)*v - 48.0*u*v*w - 18.0*v*pow(w,2.0) - 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
    sf(3,3) = -2.0*u*w - 2.0*pow(u,2.0) + v*w + pow(v,2.0);
    sf(3,4) = -2.0*pow(u,2.0) + pow(w,2.0) - 2.0*u*v + v*w;
    sf(3,5) = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf(3,6) = 1.0/2.0*(-2.0*pow(u,2.0) + pow(w,2.0) + 4.0*v*w + pow(v,2.0));
    sf(4,0) = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0) + pow(w,4.0) + 2.0*pow(u,3.0)*v + 6.0*pow(u,2.0)*v*w + 6.0*u*v*pow(w,2.0) + 2.0*v*pow(w,3.0));
    sf(4,1) = 1.0/12.0*(-2.0*pow(u,3.0) - 12.0*pow(u,2.0)*w - 18.0*u*pow(w,2.0) - 4.0*pow(w,3.0) - 6.0*pow(u,2.0)*v - 12.0*u*v*w - 6.0*v*pow(w,2.0));
    sf(4,2) = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0));
    sf(4,3) = u*w + pow(w,2.0) + u*v + v*w;
    sf(4,4) = -2.0*u*w;
    sf(4,5) = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf(4,6) = 1.0/2.0*(-pow(u,2.0) - 2.0*u*w + pow(w,2.0));
    sf(5,0) = 1.0/12.0*(2.0*u*pow(v,3.0) + pow(v,4.0)); 
    sf(5,1) = 1.0/12.0*(6.0*u*pow(v,2.0) + 2.0*pow(v,3.0)); 
    sf(5,2) = -1.0/6.0*pow(v,3.0);
    sf(5,3) = u*v; 
    sf(5,4) = 0.0;
    sf(5,5) = -1.0/2.0*pow(v,2.0); 
    sf(5,6) = -1.0/2.0*pow(v,2.0);
    sf(6,0) = 1.0/12.0*(pow(u,4.0) + 6.0*pow(u,3.0)*w + 12.0*pow(u,2.0)*pow(w,2.0) + 6.0*u*pow(w,3.0)+ pow(w,4.0) + 8.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 36.0*u*v*pow(w,2.0) + 8.0*v*pow(w,3.0) + 24.0*pow(u,2.0)*pow(v,2.0) + 60.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 24.0*u*pow(v,3.0) + 24.0*pow(v,3.0)*w + 6.0*pow(v,4.0));
   sf(6,1) = 1.0/12.0*(4.0*pow(u,3.0) + 18.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 2.0*pow(w,3.0) + 24.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 24.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w);
   sf(6,2) = 1.0/12.0*(2.0*pow(u,3.0) + 6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0) - 2.0*pow(w,3.0) + 12.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w);
   sf(6,3) = pow(u,2.0) + u*w - 2.0*v*w - 2.0*pow(v,2.0);
   sf(6,4) = -2.0*u*w - 2.0*u*v - 2.0*v*w - 2.0*pow(v,2.0);
   sf(6,5) = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
   sf(6,6) = 1.0/2.0*(pow(u,2.0) - 2.0*u*w - pow(w,2.0) - 4.0*v*w - 2.0*pow(v,2.0));
   sf(7,0) = 1.0/12.0*(pow(u,4.0) + 8.0*pow(u,3.0)*w + 24.0*pow(u,2.0)*pow(w,2.0) + 24.0*u*pow(w,3.0) + 6.0*pow(w,4.0) + 6.0*pow(u,3.0)*v + 36.0*pow(u,2.0)*v*w + 60.0*u*v*pow(w,2.0) + 24.0*v*pow(w,3.0) + 12.0*pow(u,2.0)*pow(v,2.0) + 36.0*u*pow(v,2.0)*w + 24.0*pow(v,2.0)*pow(w,2.0) + 6.0*u*pow(v,3.0) + 8.0*pow(v,3.0)*w + pow(v,4.0));
   sf(7,1) = 1.0/12.0*(2.0*pow(u,3.0) + 12.0*pow(u,2.0)*w + 12.0*u*pow(w,2.0) + 6.0*pow(u,2.0)*v - 12.0*v*pow(w,2.0) - 6.0*u*pow(v,2.0) - 12.0*pow(v,2.0)*w - 2.0*pow(v,3.0));
   sf(7,2) = 1.0/12.0*(4.0*pow(u,3.0) + 24.0*pow(u,2.0)*w + 24.0*u*pow(w,2.0) + 18.0*pow(u,2.0)*v + 48.0*u*v*w + 12.0*v*pow(w,2.0) + 12.0*u*pow(v,2.0) + 12.0*pow(v,2.0)*w + 2.0*pow(v,3.0));
   sf(7,3) = -2.0*u*w - 2.0*pow(w,2.0) - 2.0*u*v - 2.0*v*w;
   sf(7,4) = pow(u,2.0) - 2.0*pow(w,2.0) + u*v - 2.0*v*w;
   sf(7,5) = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
   sf(7,6) = 1.0/2.0*(pow(u,2.0) - 2.0*pow(w,2.0) - 2.0*u*v - 4.0*v*w - pow(v,2.0));
   sf(8,0) = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0)); 
   sf(8,1) = -1.0/6.0*pow(w,3.0); 
   sf(8,2) = 1.0/12.0*(6.0*u*pow(w,2.0) + 2.0*pow(w,3.0));
   sf(8,3) = 0.0; 
   sf(8,4) = u*w;
   sf(8,5) = -1.0/2.0*pow(w,2.0); 
   sf(8,6) = -1.0/2.0*pow(w,2.0);
   sf(9,0) = 1.0/12.0*(2.0*pow(v,3.0)*w + pow(v,4.0));
   sf(9,1) = 1.0/12.0*(6.0*pow(v,2.0)*w + 4.0*pow(v,3.0)); 
   sf(9,2) = 1.0/6.0*pow(v,3.0);
   sf(9,3) = v*w + pow(v,2.0); 
   sf(9,4) = 0.0;
   sf(9,5) = 1.0/2.0*pow(v,2.0); 
   sf(9,6) = 1.0/2.0*pow(v,2.0);
   sf(10,0) = 1.0/12.0*(2.0*u*pow(w,3.0) + pow(w,4.0) + 6.0*u*v*pow(w,2.0) + 6.0*v*pow(w,3.0) + 6.0*u*pow(v,2.0)*w + 12.0*pow(v,2.0)*pow(w,2.0) + 2.0*u*pow(v,3.0) + 6.0*pow(v,3.0)*w + pow(v,4.0));
   sf(10,1) = 1.0/12.0*(4.0*pow(w,3.0) + 18.0*v*pow(w,2.0) + 6.0*u*pow(w,2.0) + 12.0*pow(v,2.0)*w + 12.0*u*v*w + 2.0*pow(v,3.0) + 6.0*u*pow(v,2.0));
   sf(10,2) = 1.0/12.0*(2.0*pow(w,3.0) + 6.0*u*pow(w,2.0) + 12.0*v*pow(w,2.0) + 12.0*u*v*w + 18.0*pow(v,2.0)*w + 6.0*u*pow(v,2.0) + 4.0*pow(v,3.0));
   sf(10,3) = pow(w,2.0) + v*w + u*w + u*v;
   sf(10,4) = u*w + v*w + u*v + pow(v,2.0);
   sf(10,5) = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
   sf(10,6) = 1.0/2.0*(pow(w,2.0) + 4.0*v*w + 2.0*u*w + pow(v,2.0) + 2.0*u*v);
   sf(11,0) = 1.0/12.0*(pow(w,4.0) + 2.0*v*pow(w,3.0)); 
   sf(11,1) = 1.0/6.0*pow(w,3.0); 
   sf(11,2) = 1.0/12.0*(4.0*pow(w,3.0) + 6.0*v*pow(w,2.0));
   sf(11,3) = 0.0; 
   sf(11,4) = pow(w,2.0) + v*w;
   sf(11,5) = 1.0/2.0*pow(w,2.0); 
   sf(11,6) = 1.0/2.0*pow(w,2.0);
}

// for irregular patch, more subdivision is needed. 
// for different sub-element, different new nodes are selected, select-matrix (SM)
// vertex here is the original 11 vertice, so vertex is 11*3 matrix
// M(17,11); mat M1(12,17); mat M2(12,17); mat M3(12,17); mat M4(11,17);
void subdivision_matrix(mat& M, mat& SM1, mat& SM2, mat& SM3, mat& SM4){
    int N=6; double w=3.0/8.0/N; // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
    int N1=5; double w1=3.0/8.0/N1; // w1=1/N1*(5/8-(3/8+1/4*cos(2*pi/N1))^2);
    double a=3.0/8.0; double b=1.0/8.0;
    M << a << b << a << b << 0 << 0 << 0 << 0 << 0 << 0 << 0 << endr
      << b << a << a << 0 << 0 << b << 0 << 0 << 0 << 0 << 0 << endr
      << w1 << w1 << 1.0-N1*w1 << w1 << 0 << w1 << w1 << 0 << 0 << 0 << 0 << endr
      << b << 0 << a << a << 0 << 0 << b << 0 << 0 << 0 << 0 << endr
      << 0 << a << b << 0 << b << a << 0 << 0 << 0 << 0 << 0 << endr
      << 0 << b << a << 0 << 0 << a << b << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << a << b << 0 << b << a << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << b << a << 0 << 0 << a << b << 0 << 0 << 0 << endr
      << 0 << b << 0 << 0 << a << a << 0 << 0 << b << 0 << 0 << endr
      << 0 << w << w << 0 << w << 1.0-N*w <<  w << 0 << w << w << 0 << endr
      << 0 << 0 << b << 0 << 0 << a << a << 0 << 0 << b << 0 << endr
      << 0 << 0 << w << w << 0 << w << 1.0-N*w << w << 0 << w << w << endr
      << 0 << 0 << 0 << b << 0 << 0 << a << a << 0 << 0 << b << endr
      << 0 << 0 << 0 << 0 << b << a << 0 << 0 << a << b << 0 << endr
      << 0 << 0 << 0 << 0 << 0 << a << b << 0 << b << a << 0 << endr
      << 0 << 0 << 0 << 0 << 0 << b << a << 0 << 0 << a << b << endr
      << 0 << 0 << 0 << 0 << 0 << 0 << a << b << 0 << b << a << endr;
    // SM1=zeros(12,17);
    Row<int> element1; element1 << 2 << 3 << 5 << 6 << 7 << 9 << 10 << 11 << 12 << 14 << 15 << 16 << endr;
    for (int i=0; i<12; i++){
        SM1(i,element1(i)) = 1.0;
    }
    // SM2=zeros(12,17);
    Row<int> element2; element2 << 4 << 1 << 9 << 5 << 2 << 14 << 10 << 6 << 3 << 15 << 11 << 7 << endr;
    for (int i = 0; i < 12; i++){
        SM2(i,element2(i)) = 1.0;
    }
    // SM3=zeros(12,17);
    Row<int> element3; element3 << 1 << 2 << 4 << 5 << 6 << 8 << 9 << 10 << 11 << 13 << 14 << 15 << endr;
    for (int i = 0; i < 12; i++){ 
        SM3(i,element3(i)) = 1.0;
    }
    // SM4=zeros(11,17);
    Row<int> element4; element4 << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 9 << 10 << 11 << endr;
    for (int i = 0; i < 11; i++){ 
        SM4(i,element4(i)) = 1.0;
    }
}
void trans_time(mat a, mat b, mat& out){
    int rowa = a.n_rows; int cola = a.n_cols;
    int rowb = b.n_rows; int colb = b.n_cols;
    //mat out(cola,colb);
    if (rowa != rowb){
        cout<<"wrong. matrix dimentions are not correct during trans_time!"<<endl;
        exit(1);
    }else{
        for (int i = 0; i < cola; i++){
            for (int j = 0; j < colb; j++){
                double temp = 0.0;
                for (int k = 0; k < rowa; k++){
                    temp = temp + a(k,i)*b(k,j);
                }
                out(i,j) = temp;
            }
        }
    }
    //return out;
}

vec force_scale(mat force){
    int num = force.n_rows;
    vec out(num);
    #pragma omp parallel for
    for (int i = 0; i < num; i++){
        out(i) = norm(force.row(i),2);
    }
    return out;
}

rowvec tovector(mat vertex){
    int nodenum = vertex.n_rows;
    rowvec vertex_row(3*nodenum);
    #pragma omp parallel for
    for (int i = 0; i < nodenum; i++){
        vertex_row(3*i+0) = vertex(i,0);
        vertex_row(3*i+1) = vertex(i,1);
        vertex_row(3*i+2) = vertex(i,2);
    }
    return vertex_row;
}

mat tomatrix(rowvec vertex_row){
    int nodenum = vertex_row.n_cols/3;
    mat vertex(nodenum,3);
    #pragma omp parallel for
    for (int i = 0; i < nodenum; i++){
        vertex(i,0) = vertex_row(3*i+0);
        vertex(i,1) = vertex_row(3*i+1);
        vertex(i,2) = vertex_row(3*i+2);
    }
    return vertex;
}

mat setVMU(int n){
    int dotsnumber = 1;
    if (n==1){
        dotsnumber = 1;
        mat vmu(dotsnumber,3);
        vmu << 1.0/3.0 << 1.0/3.0 << 1.0/3.0 << endr;
        return vmu;
    }else if(n==2){
        dotsnumber = 3;
        mat vmu(dotsnumber,3);
        vmu << 1.0/6.0 << 1.0/6.0 << 4.0/6.0 << endr
            << 1.0/6.0 << 4.0/6.0 << 1.0/6.0 << endr
            << 4.0/6.0 << 1.0/6.0 << 1.0/6.0 << endr;
        return vmu;
    }else if(n==3){
        dotsnumber = 4;
        mat vmu(dotsnumber,3);
        vmu << 1.0/3.0 << 1.0/3.0 << 1.0/3.0 << endr
            << 1.0/5.0 << 1.0/5.0 << 3.0/5.0 << endr
            << 1.0/5.0 << 3.0/5.0 << 1.0/5.0 << endr
            << 3.0/5.0 << 1.0/5.0 << 1.0/5.0 << endr;
        /*
        vmu << 1.0/3.0 << 1.0/3.0 << 1.0/3.0 << endr
            << 2.0/15.0 << 11.0/15.0 << 2.0/15.0 << endr
            << 2.0/15.0 << 2.0/15.0 << 11.0/15.0 << endr
            << 11.0/15.0 << 2.0/15.0 << 2.0/15.0 << endr;
        */
        return vmu;
    }else if(n==4){
        dotsnumber = 6;
        mat vmu(dotsnumber,3);
        vmu << 0.44594849091597 << 0.44594849091597 << 0.10810301816807 << endr
            << 0.44594849091597 << 0.10810301816807 << 0.44594849091597 << endr
            << 0.10810301816807 << 0.44594849091597 << 0.44594849091597 << endr
            << 0.09157621350977 << 0.09157621350977 << 0.81684757298046 << endr
            << 0.09157621350977 << 0.81684757298046 << 0.09157621350977 << endr
            << 0.81684757298046 << 0.09157621350977 << 0.09157621350977 << endr;
        return vmu;
    }else if(n==5){
        dotsnumber = 7;
        mat vmu(dotsnumber,3);
        vmu << 0.33333333333333 << 0.33333333333333 << 0.33333333333333 << endr
            << 0.47014206410511 << 0.47014206410511 << 0.05971587178977 << endr
            << 0.47014206410511 << 0.05971587178977 << 0.47014206410511 << endr
            << 0.05971587178977 << 0.47014206410511 << 0.47014206410511 << endr
            << 0.10128650732346 << 0.10128650732346 << 0.79742698535309 << endr
            << 0.10128650732346 << 0.79742698535309 << 0.10128650732346 << endr
            << 0.79742698535309 << 0.10128650732346 << 0.10128650732346 << endr;
        return vmu;
    }else if(n==6){
        dotsnumber = 12;
        mat vmu(dotsnumber,3);
        vmu << 0.24928674517091 << 0.24928674517091 << 0.50142650965818 << endr
            << 0.24928674517091 << 0.50142650965818 << 0.24928674517091 << endr
            << 0.50142650965818 << 0.24928674517091 << 0.24928674517091 << endr
            << 0.06308901449150 << 0.06308901449150 << 0.87382197101700 << endr
            << 0.06308901449150 << 0.87382197101700 << 0.06308901449150 << endr
            << 0.87382197101700 << 0.06308901449150 << 0.06308901449150 << endr
            << 0.31035245103378 << 0.63650249912140 << 0.05314504984482 << endr
            << 0.63650249912140 << 0.05314504984482 << 0.31035245103378 << endr
            << 0.05314504984482 << 0.31035245103378 << 0.63650249912140 << endr
            << 0.63650249912140 << 0.31035245103378 << 0.05314504984482 << endr
            << 0.31035245103378 << 0.05314504984482 << 0.63650249912140 << endr
            << 0.05314504984482 << 0.63650249912140 << 0.31035245103378 << endr;
        return vmu;
    }
}

rowvec setVMUcoefficient(int n){
    int dotsnumber = 1;
    if (n==1){
        dotsnumber = 1;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << 1.0 << endr;
        return vmucoeff;
    }else if(n==2){
        dotsnumber = 3;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << 1.0/3.0 << 1.0/3.0 << 1.0/3.0 << endr;
        return vmucoeff;
    }else if(n==3){
        dotsnumber = 4;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << -0.56250000000000 << 0.52083333333333 << 0.52083333333333 << 0.52083333333333 << endr;
        return vmucoeff;
    }else if(n==4){
        dotsnumber = 6;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << 0.22338158967801 << 0.22338158967801 << 0.22338158967801 << 0.10995174365532 << 0.10995174365532 << 0.10995174365532 << endr;
        return vmucoeff;
    }else if(n==5){
        dotsnumber = 7;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << 0.22500000000000 << 0.13239415278851 << 0.13239415278851 << 0.13239415278851 << 0.12593918054483 << 0.12593918054483 << 0.12593918054483 << endr;
        return vmucoeff;
    }else if(n==6){
        dotsnumber = 12;
        rowvec vmucoeff(dotsnumber);
        vmucoeff << 0.11678627572638 << 0.11678627572638 << 0.11678627572638 << 0.05084490637021 << 0.05084490637021 << 0.05084490637021 << 0.08285107561837 << 0.08285107561837 << 0.08285107561837 << 0.08285107561837 << 0.08285107561837 << 0.08285107561837 << endr;
        return vmucoeff;
    }
}

void membrane_area(mat vertex, Mat<int> face, Row<int> isBoundaryFace, Mat<int> A, int GaussQuadratureN, rowvec& element_area, rowvec gqcoeff, cube shape_functions){ // A is face_ring_vertex 
    element_area.fill(0);
    
    mat VWU = setVMU(GaussQuadratureN); 
    //rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    #pragma omp parallel for
    for (int i = 0; i < face.n_rows; i++){
        if ( isBoundaryFace(i) == 1 ) 
            continue;
        
        double area = 0.0;
        {
            mat dots(A.n_cols,3); //dots(12,3); 12 nodes
            for (int j = 0; j < A.n_cols; j++){
                int nodenum = A(i,j);
                dots.row(j) = vertex.row(nodenum);
            }
            // Gaussian quadrature, 3 points
            for (int j = 0; j < VWU.n_rows; j++){
                //rowvec vwu = VWU.row(j); mat sf(12,7); shapefunctions(vwu, sf);
                mat sf = shape_functions.slice(j);
                rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                double  sqa = norm(cross(a_1,a_2),2); 
                double s = sqa; 
                area = area + 1.0/2.0*gqcoeff(j)*s; 
            }
        }
        element_area(i) = area;
    }
}

rowvec determine_spontaneous_curvature(Param param, Mat<int> face, mat vertex){
    double C0 = param.C0;
    double c0 = param.c0;
    double sigma = param.sigma; // 2*sigma is the radial region of non-zero spontaneous curvature
    bool isAdditiveScheme = param.isAdditiveScheme;
    rowvec spontcurv(face.n_rows); spontcurv.fill(c0);
    Mat<int> insertionpatch = param.insertionpatch;

    for (int i = 0; i < insertionpatch.n_rows; i++){
        for (int j = 0; j < insertionpatch.n_cols; j++){
            int facenumber = insertionpatch(i,j);
            spontcurv(facenumber) = C0;
        }
    }

    return spontcurv;
}

Row<int> determine_isinsertionpatch(Mat<int> face, Mat<int> insertionpatch){
    int facenumber = face.n_rows;
    Row<int> out(facenumber);
    for (int i = 0; i < facenumber; i++){
        bool Isinsertionpatch = false;
        for (int j = 0; j < insertionpatch.n_rows; j++){
            for (int m = 0; m < insertionpatch.n_cols; m++){
                if ( i == insertionpatch(j,m) )
                   Isinsertionpatch = true; 
            }
        }
        if ( Isinsertionpatch == true ){
            out(i) = 1;
        }else{
            out(i) = 0;
        }
    }
    return out;
}

mat manage_ghost_force(mat force, Mat<int> face, Row<int> isGhostFace, Row<int> isGhostVertex, Param param){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3)/2 * a; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;
    mat fxyz = force;
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( isGhostFace(i) != 1 )
                continue;
            int node1 = face(i,0);
            int node2 = face(i,1);
            int node3 = face(i,2);
            fxyz.row(node1) = force.row(node1) * 0.0;
            fxyz.row(node2) = force.row(node2) * 0.0;
            fxyz.row(node3) = force.row(node3) * 0.0;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        #pragma omp parallel for 
        for (int i = 0; i < vertexnum; i++){
            if ( isGhostVertex(i) == 1 )
                fxyz.row(i) = force.row(i) * 0.0;
        }
    }else if ( param.isBoundaryFree == true ){   
        for ( int i = 0; i < n+1; i++ ){
            int j = 0;
            int index = (n+1)*j + i;
            fxyz.row(index) = force.row(index) * 0.0;
            index = (n+1)*m + i;
            fxyz.row(index) = force.row(index) * 0.0;
        }
        for ( int j = 0; j < m+1; j++ ){
            int i = 0;
            int index = (n+1)*j + i;
            fxyz.row(index) = force.row(index) * 0.0;
            index = (n+1)*j + n;
            fxyz.row(index) = force.row(index) * 0.0;
        }
    }
    return fxyz;
}

mat update_vertex(mat vertex, Mat<int> face, double a, mat force, Param param, Row<int> isGhostFace, Row<int> isGhostVertex){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double aa = dx;
    double dy = sqrt(3)/2 * aa; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;
    
    mat vertexnew = vertex + a*force;
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( isGhostFace(i) != 1 )
                continue;
            int node1 = face(i,0);
            int node2 = face(i,1);
            int node3 = face(i,2);
            vertexnew.row(node1) = vertex.row(node1);
            vertexnew.row(node2) = vertex.row(node2);
            vertexnew.row(node3) = vertex.row(node3);
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
            vertexnew.row(index0) = vertex.row(index0) + ( vertexnew.row(index00)-vertex.row(index00) );
            vertexnew.row(index1) = vertex.row(index1) + ( vertexnew.row(index11)-vertex.row(index11) );
            vertexnew.row(index2) = vertex.row(index2) + ( vertexnew.row(index22)-vertex.row(index22) );
        }
        #pragma omp parallel for 
       for (int i = 0; i < n; i++){
            int index0 = n*(m-2) + i;
            int index1 = n*(m-1) + i;
            int index2 = n*(m-0) + i;
            int index00 = n*4 + i;
            int index11 = n*5 + i;
            int index22 = n*6 + i;
            vertexnew.row(index0) = vertex.row(index0) + ( vertexnew.row(index00)-vertex.row(index00) );
            vertexnew.row(index1) = vertex.row(index1) + ( vertexnew.row(index11)-vertex.row(index11) );
            vertexnew.row(index2) = vertex.row(index2) + ( vertexnew.row(index22)-vertex.row(index22) );
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
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        }
        index1 = (n+1)*1+0; 
        index2 = (n+1)*2+0;
        index3 = (n+1)*1+1;
        index4 = (n+1)*2+1;
        vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        // right side
        index1 = (n+1)*1 + n; 
        index2 = (n+1)*2 + n-1;
        index3 = (n+1)*1 + n-1;
        index4 = (n+1)*2 + n-2;
        vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
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
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        }
        // bottom 
        for (int i = 0; i < n; i++){
            index1 = (n+1)*0 + i;
            index2 = (n+1)*1 + i;
            index3 = (n+1)*1 + i+1;
            index4 = (n+1)*2 + i;
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        }
        index1 = (n+1)*0 + n; 
        index2 = (n+1)*0 + n-1;
        index3 = (n+1)*1 + n;
        index4 = (n+1)*1 + n-1;
        vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        // top 
        for (int i = 1; i < n+1; i++){
            index1 = (n+1)*m + i;
            index2 = (n+1)*(m-1) + i;
            index3 = (n+1)*(m-1) + i+1;
            index4 = (n+1)*(m-2) + i;
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
        }
        index1 = (n+1)*m + n; 
        index2 = (n+1)*m + n-1;
        index3 = (n+1)*(m-1) + n;
        index4 = (n+1)*(m-1) + n-1;
        vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
    }
    
    return vertexnew;
}

