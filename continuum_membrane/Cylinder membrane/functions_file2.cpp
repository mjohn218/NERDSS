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
    bool   isBoundaryFixed = false;
    bool   isBoundaryPeriodic = false;
    bool   isBoundaryFree = false;
    bool   isBoundaryForced = false;
    double extforce = 0.0;
    double Radius;
    double Length;
    double l;     // for flat subdivision 
}; 

void finddot(int a, int b, Mat<int> vertexi_face, Mat<int> face, Row<int>& dot){
    Row<int> aface; aface = vertexi_face.row(a);
    Row<int> bface; bface = vertexi_face.row(b);
    Row<int> A(2); // A has two element, one of which is face_i
    int shu= -1;
    for (int i=0; i<aface.n_cols; i++){
        for (int j=0; j<bface.n_cols; j++){
            if (aface(i) == bface(j) && aface(i) != -1){
                shu = shu + 1;
                A(shu) = aface(i);
                if (aface(i) == -1){
                    cout<<"wrong to finddot!!"<<endl;
                    exit(1);
                }
            }
        }
    }
    for (int j=0; j<2; j++){
        for (int k=0; k<3; k++) {
            if (face(A(j),k)!=a && face(A(j),k)!=b){ // get out the vertex in face(A(j)) that is not a nor b
                dot(j) = face(A(j),k);
            }
        }
    }
}

Row<int> find_two_faces(int node1, int node2, Mat<int> vertexi_face){
    Row<int> twofaces(2); twofaces.fill(-1);
    int number = 0;
    for (int i = 0; i < 6; i ++){
        if ( vertexi_face(node1,i) == -1 ) continue;
        for (int j = 0; j < 6; j++){
            if ( vertexi_face(node2,j) == -1 ) continue;
            if ( vertexi_face(node1,i) == vertexi_face(node2,j) ){
                twofaces(number) = vertexi_face(node1,i);
                number++;
            }
        }
    }
    if (twofaces(0) == -1 || twofaces(1) == -1){
        cout<<"Wrong: in find_two_faces, unsuccessful! "<<endl;
        exit(0);
    }
    return twofaces;
}


bool isSame(mat A, mat B){
    bool out = true;
    for (int i=0; i<A.n_rows; i++){
        for (int j=0; j<A.n_cols; j++){
            if (A(i,j) != B(i,j))
                out = false;
        }
    }
    return out;
}

Mat<double> setvertex_Loop_scheme(Param& param){ // vertex position
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;

    int n = round(2.0*M_PI*R/l); 
    R = l/2.0/sin(M_PI/n); //param.Radius = R;
    cout << "Cylinder radius = "<< R <<" nm"<< endl;
    double dtheta = 2.0*M_PI/n; // circle division
    double a = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * a; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    L = m*dz; //param.Length = L;

    if (param.isBoundaryFixed == true){
        cout << "FixedBoundary Cylinder length = "<< L-2.0*dz <<" nm"<< endl;
    }else if (param.isBoundaryPeriodic == true){
        cout << "PeriodicBoundary Cylinder length = "<< L-6.0*dz <<" nm"<< endl;
    }
    
    int nodenum = n*(m+1);
    mat vertex(nodenum,3); vertex.fill(0.0);
    #pragma omp parallel for
    for (int j = 0; j < m+1; j++){
        bool isEvenJ = false;
        if ( pow(-1.0,j) > 0.0 ){
            isEvenJ = true;
        }
        double z = j*dz;
        for (int i = 0; i < n; i++){
            int index = n*j + i;
            double theta = i*dtheta;
            if ( isEvenJ == false ){
                theta = theta + dtheta/2.0;
            }
            double x = R*cos(theta);
            double y = R*sin(theta);
            vertex(index,0) = x; 
            vertex(index,1) = y;
            vertex(index,2) = z;
        }
    }
    return vertex;
}
Mat<int> setface_Loop_scheme(Param param){ // face and its surrounding vertex
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;

    int n = round(2.0*M_PI*R/l); 
    R = l/2.0/sin(M_PI/n);
    double dtheta = 2.0*M_PI/n; // circle division
    double a = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * a; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    
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
            if ( isEvenJ == true ){
                node1 = n*j + i;
                node2 = n*(j+1) + i;
                node3 = n*j + (i+1);
                node4 = n*(j+1) + (i+1);
                if ( i+1 == n ){
                    node3 = n*j + 0;
                    node4 = n*(j+1) + 0;
                }
            }else{
                node1 = n*(j+1) + i;
                node2 = n*(j+1) + (i+1);
                node3 = n*j + i;
                node4 = n*j + (i+1);
                if ( i+1 == n ){
                    node2 = n*(j+1) + 0;
                    node4 = n*j + 0;
                }
            }
            face(index,0) = node1; face(index,1) = node2; face(index,2) = node3;
            face(index+1,0) = node2; face(index+1,1) = node3; face(index+1,2) = node4;
        }
    }
    return face;
}

Row<int> isEndVertex(Param param){
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;

    int n = round(2.0*M_PI*R/l); 
    R = l/2.0/sin(M_PI/n);
    double dtheta = 2.0*M_PI/n; // circle division
    double a = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * a; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    
    int nodenum = n*(m+1);
    Row<int> isBoundaryNode(nodenum); isBoundaryNode.fill(-1);
    int j = 0;
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        int index = n*j + i;
        isBoundaryNode(index) = 1; // if element is 1, then this vertex is on boundary
    }
    j = m;
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
         int index = n*j + i;
        isBoundaryNode(index) = 1; // if element is 1, then this vertex is on boundary
    }
    return isBoundaryNode;
}

Row<int> isBoundaryface(Mat<int> face, mat vertex, Param param){
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;

    int n = round(2.0*M_PI*R/l); 
    R = l/2.0/sin(M_PI/n);
    double dtheta = 2.0*M_PI/n; // circle division
    double a = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * a; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    L = m*dz;

    int facenum = m*n*2;
    Row<int> isBoundaryFace(facenum); isBoundaryFace.fill(-1);
    double ztop, zbottom;
    if (param.isBoundaryFixed == true){
        ztop = L - dz;
        zbottom = 0.0 + dz;
    }else if (param.isBoundaryPeriodic == true){
        ztop = L - 3.0*dz;
        zbottom = 0.0 + 3.0*dz;
    }
    #pragma omp parallel for 
    for (int i = 0; i < facenum; i++){
        int node1 = face(i,0);
        int node2 = face(i,1);
        int node3 = face(i,2);
        double z = 1.0/3.0 * ( vertex(node1,2) + vertex(node2,2) + vertex(node3,2) );
        if ( z < zbottom || z > ztop ){
            isBoundaryFace(i) = 1; // if element is not -1, then this face is on boundary
        }
        if ( z < dz || z > L - dz ){ // the end faces
            isBoundaryFace(i) = 2; // if element is 2, then this face is on the end
        }
    }
    return isBoundaryFace;
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

// find out the one-ring vertices aound face_i. It could be 12 or less. 
Mat<int> one_ring_vertices(Mat<int> face, Mat<int> face_with_vertexi, Row<int> isBoundaryFace){ 
    int facenum = face.n_rows;
    Mat<int>  ring_vertices(facenum,12);  ring_vertices.fill(-1);
    #pragma omp parallel for
    for (int i = 0; i < facenum; i++){
        if ( isBoundaryFace(i) == 2 ){ // exclude those end faces
            continue;
        }
        int d1, d2, d3, d5, d6, d9, d10, d11, d12;
        int d4 = face(i,0); 
        int d7 = face(i,1);
        int d8 = face(i,2);
        Row<int> dot3(2); finddot(d4,d7,face_with_vertexi,face,dot3);
        Row<int> dot5(2); finddot(d4,d8,face_with_vertexi,face,dot5);
        Row<int> dot11(2); finddot(d7,d8,face_with_vertexi,face,dot11);
        for (int j=0; j<2; j++){
            if ( dot3(j) != d8 )  d3 = dot3(j);
            if ( dot5(j) != d7 )  d5 = dot5(j);
            if ( dot11(j) != d4 ) d11 = dot11(j);
        }
        Row<int> dot1(2); finddot(d4,d3,face_with_vertexi,face,dot1);
        Row<int> dot2(2); finddot(d4,d5,face_with_vertexi,face,dot2);
        Row<int> dot6(2); finddot(d7,d3,face_with_vertexi,face,dot6);
        Row<int> dot10(2); finddot(d7,d11,face_with_vertexi,face,dot10);
        Row<int> dot9(2); finddot(d8,d5,face_with_vertexi,face,dot9);
        Row<int> dot12(2); finddot(d8,d11,face_with_vertexi,face,dot12);
        for (int j = 0; j < 2; j++){
            if ( dot1(j) != d7 ) d1 = dot1(j);  
            if ( dot2(j) != d8 ) d2 = dot2(j); 
            if ( dot6(j) != d4 ) d6 = dot6(j); 
            if ( dot10(j) != d8 ) d10 = dot10(j);  
            if ( dot9(j) != d4 ) d9 = dot9(j); 
            if ( dot12(j) != d7 ) d12 = dot12(j); 
        }
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
        if ( isBoundaryFace(i) > -1 ) // exclude those boundary and end faces
            continue;
        
        double area = 0.0;
        //if  ( A(i,0) != A(i,1) )
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

rowvec determine_spontaneous_curvature(Param param, Mat<int> face, mat vertex, Mat<int> insertionpatch){
    double C0 = param.C0;
    double c0 = param.c0;
    double sigma = param.sigma; // 2*sigma is the radial region of non-zero spontaneous curvature
    bool isAdditiveScheme = param.isAdditiveScheme;
    rowvec spontcurv(face.n_rows); spontcurv.fill(c0);
    
    for (int i = 0; i < insertionpatch.n_rows; i++){
        for (int j = 0; j < insertionpatch.n_cols; j++){
            int facenumber = insertionpatch(i,j);
            spontcurv(facenumber) = C0;
        }
    }
    if ( sigma < 1e-9 ) 
        return spontcurv;
    
    if (insertionpatch.n_rows < 2){
        for (int i = 0; i < face.n_rows; i++ ){
            int j = 0;
            rowvec centeri = 1.0/3.0 * ( vertex.row(face(i,0)) + vertex.row(face(i,1)) + vertex.row(face(i,2)) );
            double dismin = 2.0*param.Radius;
            for (int k = 0; k < insertionpatch.n_cols; k++){
                int facenumber = insertionpatch(j,k);
                rowvec center0 = 1.0/3.0 * ( vertex.row(face(facenumber,0)) + vertex.row(face(facenumber,1)) + vertex.row(face(facenumber,2)) );
                rowvec disvec = centeri - center0; 
                double distance = norm(disvec,2.0);
                if (distance < dismin){
                    dismin = distance;
                }
            }
            spontcurv(i) = - abs(C0) * exp( - pow(dismin/sigma,2.0)/2.0 ); // here the spont_curvature is negative sign
            if ( abs(spontcurv(i)) < 1.0e-15 )
                spontcurv(i) = 0.0;
        }
        return spontcurv;
    }
    
    for (int i = 0; i < face.n_rows; i++ ){
        rowvec centeri = 1.0/3.0 * ( vertex.row(face(i,0)) + vertex.row(face(i,1)) + vertex.row(face(i,2)) );
        rowvec spontcurvs(insertionpatch.n_rows); spontcurvs.fill(0.0);
        for (int j = 0; j < insertionpatch.n_rows; j++){
            double dismin = 2.0*param.Radius;
            for (int k = 0; k < insertionpatch.n_cols; k++){
                int facenumber = insertionpatch(j,k);
                rowvec center0 = 1.0/3.0 * ( vertex.row(face(facenumber,0)) + vertex.row(face(facenumber,1)) + vertex.row(face(facenumber,2)) );
                rowvec disvec = centeri - center0;
                double distance = norm(disvec,2.0);
                if (distance < dismin){
                    dismin = distance;
                }
            }
            spontcurvs(j) = abs(C0) * exp( - pow(dismin/sigma,2.0)/2.0 );
        }
        if (isAdditiveScheme == true){
            /*
            double shu = 0.0;
            double H = 1.0/R;
            for (int j = 0; j < insertionpatch.n_rows; j++){
                shu = shu + pow(2.0*H-spontcurvs(j),2.0) - pow(2.0*H-0.0,2.0); // sum of the membrane energy change by all the insertions
            }
            shu = shu +  pow(2.0*H-0.0,2.0);
            if (shu < 0.0){ 
                cout<<"Note: when decide the enhanced spontaneous curvature, nonlinear effect happens. Then no additive scheme is utilized! "<<endl;
                spontcurv(i) = - max(spontcurvs);
            }else{
                if (max(spontcurvs) >= 2.0*H){
                    spontcurv(i) = -(2.0*H + sqrt(shu)); // here the spont_curvature is negative sign
                }else{
                    spontcurv(i) = -(2.0*H - sqrt(shu)); 
                }
                if ( spontcurv(i) > 0.0 ){
                    cout<<"Wrong: not efficient spontaneous curvature is calculated! Exit!"<<endl;
                    exit(0);
                }
            }
             */
        }else{
            spontcurv(i) = - max(spontcurvs); // here the spont_curvature is negative sign
        }
        
        if ( abs(spontcurv(i)) < 1.0e-15 )
            spontcurv(i) = 0.0;
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

mat manage_ghost_force(mat force, Mat<int> face, Row<int> isBoundaryFace, Param param){
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;
    int n = round(2.0*M_PI*R/l);
    R = l/2.0/sin(M_PI/n);
    double dtheta = 2.0*M_PI/n; // circle division
    double a = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * a; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    
    int vertexnum = n*(m+1);
    int facenum = m*n*2;
    mat fxyz = force;
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( isBoundaryFace(i) < 1 )
                continue;
            int node1 = face(i,0);
            int node2 = face(i,1);
            int node3 = face(i,2);
            fxyz.row(node1) = fxyz.row(node1) * 0.0;
            fxyz.row(node2) = fxyz.row(node2) * 0.0;
            fxyz.row(node3) = fxyz.row(node3) * 0.0;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        /*
        mat down0(n,3), down1(n,3), upM_1(n,3), upM(n,3);
        // find out the forces on boundaries, each up and down lines
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index1 = n*0 + i;
            int index2 = n*1 + i;
            down0.row(i) = force.row(index1);
            down1.row(i) = force.row(index2);
            int index3 = n*(m-1) + i;
            int index4 = n*m + i;
            upM_1.row(i) = force.row(index3);
            upM.row(i) = force.row(index4);
        }
        // add the boundary forces to the conrresponding inner vertices
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index1 = n*1 + i;
            int index2 = n*2 + i;
            fxyz.row(index1) = fxyz.row(index1) + upM_1.row(i);
            fxyz.row(index2) = fxyz.row(index2) + upM.row(i);
            int index3 = n*(m-2) + i;
            int index4 = n*(m-1) + i;
            fxyz.row(index3) = fxyz.row(index3) + down0.row(i);
            fxyz.row(index4) = fxyz.row(index4) + down1.row(i);
        }
        */
       #pragma omp parallel for
       for (int i = 0; i < n; i++){
           for (int j = m-2; j < m+1; j ++){
               int index = n*j + i;
               fxyz.row(index) = fxyz.row(index) * 0.0;
           }
           for (int j = 0; j < 3; j ++){
               int index = n*j + i;
               fxyz.row(index) = fxyz.row(index) * 0.0;
           }
       }
    }else if ( param.isBoundaryFree == true ){
          
    }else if ( param.isBoundaryForced == true ){
        #pragma omp parallel for
        for (int j = 0; j <= 1; j++){
            for (int i = 0; i < n; i++){
                int index = n*j + i;
                rowvec extforce(3); extforce << 0.0 << 0.0 << - param.extforce << endr;
                fxyz.row(index) = fxyz.row(index) + extforce;
            }
        }
        for (int j = m-1; j <= m; j++){
            for (int i = 0; i < n; i++){
                int index = n*j + i;
                rowvec extforce(3); extforce << 0.0 << 0.0 <<  param.extforce << endr;
                fxyz.row(index) = fxyz.row(index) + extforce;
            }
        }
    }

    return fxyz;
}

mat update_vertex(mat vertex, Mat<int> face, double a, mat force, Param param, Row<int> isBoundaryFace){
    double R = param.Radius;
    double L = param.Length;
    double l = param.l;
    int n = round(2.0*M_PI*R/l); 
    R = l/2.0/sin(M_PI/n);
    double dtheta = 2.0*M_PI/n; // circle division
    double aa = 2.0*R*sin(dtheta/2.0); // side of the triangle
    double dz = sqrt(3.0)/2.0 * aa; // z axis division
    int m = round(L/dz);           // z axis division
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; } // m must be an even number, easier for boundary setup
    L = m*dz; 
    if (param.isBoundaryFixed == true){
        L = m*dz - 2.0*dz;
    }else if (param.isBoundaryPeriodic == true){
        L = m*dz - 6.0*dz;
    }

    int vertexnum = n*(m+1);
    int facenum = m*n*2;
    
    mat vertexnew = vertex + a*force;
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( isBoundaryFace(i) < 1 )
                continue;
            int node1 = face(i,0);
            int node2 = face(i,1);
            int node3 = face(i,2);
            vertexnew.row(node1) = vertex.row(node1);
            vertexnew.row(node2) = vertex.row(node2);
            vertexnew.row(node3) = vertex.row(node3);
        }
    }else if ( param.isBoundaryPeriodic == true ){
        /*
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index1 = n*0 + i;
            int index2 = n*(m-2) + i;
            vertexnew.row(index1) = vertex.row(index1) + ( vertexnew.row(index2)-vertex.row(index2) );
        }
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index1 = n*m + i;
            int index2 = n*2 + i;
            vertexnew.row(index1) = vertex.row(index1) + ( vertexnew.row(index2)-vertex.row(index2) );
        }
        */
       #pragma omp parallel for 
       for (int i = 0; i < n; i++){
            int index0 = n*0 + i;
            int index1 = n*1 + i;
            int index2 = n*2 + i;
            int index00 = n*(m-6) + i;
            int index11 = n*(m-5) + i;
            int index22 = n*(m-4) + i;
            /*
            vertexnew.row(index0) = vertex.row(index00) - L;
            vertexnew.row(index1) = vertex.row(index11) - L;
            vertexnew.row(index2) = vertex.row(index22) - L;
            */
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
            /*
            vertexnew.row(index0) = vertex.row(index00) + L;
            vertexnew.row(index1) = vertex.row(index11) + L;
            vertexnew.row(index2) = vertex.row(index22) + L;
            */
            vertexnew.row(index0) = vertex.row(index0) + ( vertexnew.row(index00)-vertex.row(index00) );
            vertexnew.row(index1) = vertex.row(index1) + ( vertexnew.row(index11)-vertex.row(index11) );
            vertexnew.row(index2) = vertex.row(index2) + ( vertexnew.row(index22)-vertex.row(index22) );
        }
    }else if ( param.isBoundaryFree == true ){
        int index1, index2, index3, index4;
        for (int i = 0; i < n; i++){
            index1 = n*0 + i;
            index2 = n*1 + i-1;
            index3 = n*1 + i;
            index4 = n*2 + i;
            if ( i < 1 ){
                index2 = n*1 + n-1;
            }
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
            /*
            // move the ghost vertex onto the cylinder surface
            double x = vertexnew(index1,0);
            double y = vertexnew(index1,1);
            vertexnew(index1,0) = x * R/sqrt(x*x + y*y);
            vertexnew(index1,1) = y * R/sqrt(x*x + y*y);
            */
            index1 = n*m + i;
            index2 = n*(m-1) + i-1;
            index3 = n*(m-1) + i;
            index4 = n*(m-2) + i;
            if ( i < 1 ){
                index2 = n*(m-1) + n-1;
            }
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
            /*
            x = vertexnew(index1,0);
            y = vertexnew(index1,1);
            vertexnew(index1,0) = x * R/sqrt(x*x + y*y);
            vertexnew(index1,1) = y * R/sqrt(x*x + y*y);
            */
        }
    }else if ( param.isBoundaryForced == true ){
        int index1, index2, index3, index4;
        for (int i = 0; i < n; i++){
            index1 = n*0 + i;
            index2 = n*1 + i-1;
            index3 = n*1 + i;
            index4 = n*2 + i;
            if ( i < 1 ){
                index2 = n*1 + n-1;
            }
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
            /*
            // move the ghost vertex onto the cylinder surface
            double x = vertexnew(index1,0);
            double y = vertexnew(index1,1);
            vertexnew(index1,0) = x * R/sqrt(x*x + y*y);
            vertexnew(index1,1) = y * R/sqrt(x*x + y*y);
            */
            index1 = n*m + i;
            index2 = n*(m-1) + i-1;
            index3 = n*(m-1) + i;
            index4 = n*(m-2) + i;
            if ( i < 1 ){
                index2 = n*(m-1) + n-1;
            }
            vertexnew.row(index1) = vertexnew.row(index2) + vertexnew.row(index3) - vertexnew.row(index4);
            /*
            x = vertexnew(index1,0);
            y = vertexnew(index1,1);
            vertexnew(index1,0) = x * R/sqrt(x*x + y*y);
            vertexnew(index1,1) = y * R/sqrt(x*x + y*y);
            */
        }
    }
    
    return vertexnew;
}
