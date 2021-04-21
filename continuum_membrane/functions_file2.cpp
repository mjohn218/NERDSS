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
    double V0;
    double S; // area, total area 
    double V;
    double C0; // spontaneous curvature of insertions
    double c0; // spontaneous curvature of membrane
    double meanL;
    double gama_shape;
    double gama_area;
    double R;
    double sigma = 0.0;
    int    GaussQuadratureN;
    bool isInsertionAreaConstraint = false;
    bool isAdditiveScheme = false;
    double s0;
};

void seticosahedron(mat& vertex, Mat<int>& face, double r){
    // build the vertex
    rowvec t(5);
    for (int i=0; i<5; i++){ t(i) = i * (2.0*M_PI/5.0); }
    mat vertex1(5,3);  // vertex of up pentagon
    for (int i=0; i<5; i++){ 
        rowvec v;
        v << cos(t(i)) << sin(t(i)) << 0.0;
        vertex1.row(i) = v;
    }
    t = t + M_PI/5.0;
    double a = 2.0*sin(M_PI/5.0);   // side length
    mat vertex2(5,3);    // vertex of down pentagon
    for (int i=0; i<5; i++){ 
        rowvec v;
        v << cos(t(i)) << sin(t(i)) << -a*sqrt(3.0)/2.0;
        vertex2.row(i) = v;
    }
    double h = sqrt(a*a-1.0);          // distance between upest point and pentagon
    rowvec vertex3; vertex3<<0.0<<0.0<<h; // upest vertex
    rowvec vertex4; vertex4<<0.0<<0.0<<-a*sqrt(3.0)/2.0-h; // downest vertex
    //mat vertex(12,3);
    vertex.row(0) = vertex3;
    vertex.rows(1,5) = vertex1;
    vertex.rows(6,10) = vertex2;
    vertex.row(11) = vertex4;
    vertex.col(2) += a*sqrt(3.0)/4.0; // move the center to the origin
    for (int i=0; i<12; i++) {
        vertex.row(i) = vertex.row(i)/norm(vertex.row(i),2)*r;  // move the vertex to the surface of the sphere 
    }
    // build the face
    face<< 0 << 1 << 2 << endr
        << 0 << 2 << 3 << endr
        << 0 << 3 << 4 << endr
        << 0 << 4 << 5 << endr
        << 0 << 5 << 1 << endr
        << 1 << 6 << 2 << endr
        << 2 << 7 << 3 << endr
        << 3 << 8 << 4 << endr
        << 4 << 9 << 5 << endr
        << 5 << 10 << 1 << endr
        << 2 << 6 << 7 << endr
        << 3 << 7 << 8 << endr
        << 4 << 8 << 9 << endr
        << 5 << 9 << 10 << endr
        << 1 << 10 << 6 << endr
        << 11 << 7 << 6 << endr
        << 11 << 8 << 7 << endr
        << 11 << 9 << 8 << endr
        << 11 << 10 << 9 << endr
        << 11 << 6 << 10 << endr;
}  

// find the faces around one vertex
void vertexi_face_with_it(mat vertex, Mat<int> face, Mat<int>& nodei_face_with_it){
    #pragma omp parallel for // private(facenumber)
     for (int i=0; i<vertex.n_rows; i++){
        int facenumber = -1;
        for (int j=0; j<face.n_rows; j++){
            for (int k=0; k<face.n_cols; k++){
                if ( i == face(j,k) ){
                    facenumber = facenumber + 1;
                    nodei_face_with_it(i,facenumber)=j;
                }
            }
        }  
     }
}
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

void neighbor_vertices(Mat<int> vertexi_faces, Mat<int> face, Mat<int>& nearby_vertices){    
    #pragma omp parallel for
    for (int i = 0; i < vertexi_faces.n_rows; i++){
        //how many neighbor vertices
        int N = 6;
        if ( vertexi_faces(i,5) == -1 ){
            N = 5;
        }
        Row<int> A(N); A.fill(-1);
        int shu = -1;
        for (int j = 0; j < N; j++){
            int Numface = vertexi_faces(i,j); 
            for (int k = 0; k < 3; k++){
                if ( face(Numface,k) != i ){
                    bool islisted = false;
                    for (int m = 0; m < N; m++){
                        if (face(Numface,k) == A(m)){
                            islisted = true;
                        }
                    }
                    if ( islisted == false ){
                         shu = shu + 1;
                         A(shu) = face(Numface,k);
                    }
                }
            }
        }
        for (int j = 0; j < N; j++){
            nearby_vertices(i,j) = A(j);
        }
    }
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

void one_ring_vertices(Mat<int> face,Mat<int> face_with_vertexi, Mat<int>& ring_vertices){   
    #pragma omp parallel for
for (int i=0; i<face.n_rows; i++){
    int d1, d2, d3, d5, d6, d9, d10, d11, d12;
    int d4 = face(i,0); 
    int d7 = face(i,1);
    int d8 = face(i,2);
    // for the irregular patch, make the irregular vertex as d4
    if ( face_with_vertexi(face(i,1),5) == -1 ){
        d4 = face(i,1); d7 = face(i,2); d8 = face(i,0);
    }
    if ( face_with_vertexi(face(i,2),5) == -1 ){
        d4 = face(i,2); d7 = face(i,0); d8 = face(i,1);
    }
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
}

void setsphere_Loop_scheme(mat& vertex, Mat<int>& face, double& meanl, double r, double l){
    double a = r*2.0*sin(M_PI/5.0)/(sqrt(4.0*pow(sin(M_PI/5.0),2.0)-1.0)+0.5*cos(M_PI/5.0)); // a is the side length of icosahedron
    int n = round(log(a/l)/log(2.0));                                // division times to make the side-length as l.
    seticosahedron(vertex, face, r);                               // set an icosahedron for division into a sphere   
    for (int j = 0; j < n; j++){
        mat oldvertex; oldvertex = vertex;
        Mat<int> oldface; oldface = face;
        Mat<int> vertexi_face(vertex.n_rows,6);  vertexi_face.fill(-1); // all faces that has this vertex_i
        vertexi_face_with_it(oldvertex, oldface, vertexi_face); // if vertexi_face(x,5)=-1, means this vertex has only 5 nearby faces
        Mat<int> vertexi_nearby(vertex.n_rows, 6); vertexi_nearby.fill(-1); // all nearby vertices around this vertex_i, 5 or 6.
        neighbor_vertices(vertexi_face, oldface, vertexi_nearby);  
        int facenumber = oldface.n_rows;
        Mat<int> newface(facenumber*3, 3);
        for (int i=0; i < facenumber; i++){
            int vertexnumber = vertex.n_rows;
            ///////////////////////////////////////////////////////
            //  new vertices
            Row<int> dot1(2); finddot(oldface(i,0),oldface(i,1),vertexi_face,oldface, dot1);
            rowvec newvertex1(3); 
                   newvertex1 = 3.0/8.0*(vertex.row(oldface(i,0))+vertex.row(oldface(i,1)))+1.0/8.0*(vertex.row(dot1(0))+vertex.row(dot1(1)));
            Row<int> dot2(2); finddot(oldface(i,1),oldface(i,2),vertexi_face,oldface, dot2);
            rowvec newvertex2(3);
                   newvertex2 = 3.0/8.0*(vertex.row(oldface(i,1))+vertex.row(oldface(i,2)))+1.0/8.0*(vertex.row(dot2(0))+vertex.row(dot2(1)));
            Row<int> dot3(2); finddot(oldface(i,2),oldface(i,0),vertexi_face,oldface, dot3);
            rowvec newvertex3(3);
                   newvertex3 = 3.0/8.0*(vertex.row(oldface(i,0))+vertex.row(oldface(i,2)))+1.0/8.0*(vertex.row(dot3(0))+vertex.row(dot3(1)));
            // check whether the new vertexs are actually those existed ones
            int n1=1; int n2=1; int n3=1;
            int new1, new2, new3;
            for (int k=0; k<vertexnumber; k++) { 
                if ( isSame(newvertex1,vertex.row(k)) ){
                    n1=0;
                    new1=k;
                }
                if ( isSame(newvertex2,vertex.row(k)) ){
                    n2=0;
                    new2=k;
                }
                if ( isSame(newvertex3,vertex.row(k)) ){
                    n3=0;
                    new3=k;
                }
            }
            if (n1==1){                     // if the new vertex is totally new,not same as existed one
                new1=n1+vertexnumber-1;       // the new vertex's number 
                int k = vertex.n_rows;
                vertex.insert_rows(k, newvertex1); // add the new vertex to matrix
            }
            if (n2==1){
                new2=n1+n2+vertexnumber-1;
                int k = vertex.n_rows;
                vertex.insert_rows(k,newvertex2);
            }
            if (n3==1){
                new3=n1+n2+n3+vertexnumber-1;
                int k = vertex.n_rows;
                vertex.insert_rows(k,newvertex3);
            }
            /////////////////////////////////////////////////
            // new face
            Row<int> temp(3); 
            temp << oldface(i,0) << new1 << new3;
            newface.row(3*i) = temp ;
            temp << oldface(i,1) << new2 << new1;
            newface.row(3*i+1) = temp;
            temp << oldface(i,2) << new3 << new2;
            newface.row(3*i+2) = temp;
            temp << new1 << new2 << new3;
            face.row(i) = temp;
        }
        // update the positions of oldvertices
        for (int i = 0; i < oldvertex.n_rows; i++){
            int N = 6;
            if ( vertexi_nearby(i,5) < 0 ){
                N = 5;
            }
            // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
            double w = 3.0/8.0/N; // w=1/(N+3/8/w);
            oldvertex.row(i)=(1.0 - N*w) * oldvertex.row(i);
            for (int k = 0; k < N; k++){
                oldvertex.row(i) = oldvertex.row(i) + w * oldvertex.row(vertexi_nearby(i,k));
            }
        }
        // update the matrix for face and vertex.
        oldface = face;
        face.set_size(oldface.n_rows + newface.n_rows, oldface.n_cols);
        face.rows(0,oldface.n_rows-1) = oldface;
        face.rows(oldface.n_rows, oldface.n_rows + newface.n_rows-1) = newface;
        vertex.rows(0, oldvertex.n_rows-1) = oldvertex;
    }
    //////////////////////////////////////////////////////////
    // calculate the smallest and largest side-length; out put smallest and largest 
    vec sidelength(3*face.n_rows); 
    #pragma omp parallel for
    for (int i = 0; i < face.n_rows; i++){
        sidelength(3*i) = norm(vertex.row(face(i,0))-vertex.row(face(i,1)),2);
        sidelength(3*i+1) = norm(vertex.row(face(i,0))-vertex.row(face(i,2)),2);
        sidelength(3*i+2) = norm(vertex.row(face(i,2))-vertex.row(face(i,1)),2);
    }
    meanl = mean(sidelength);
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

void cell_area_volume(mat vertex, Mat<int> face, Mat<int> A, int GaussQuadratureN, rowvec& element_area, rowvec& element_volume, rowvec gqcoeff, cube shape_functions){ // A is face_ring_vertex 
    element_area.fill(0);
    element_volume.fill(0);
    
    mat VWU = setVMU(GaussQuadratureN); 
    //rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    #pragma omp parallel for
    for (int i = 0; i < face.n_rows; i++){
        double area = 0.0;
        double volume = 0.0;
        if  ( A(i,0) != A(i,1) ){
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
                rowvec d(3); d = cross(a_1,a_2)/sqa; 
                double s = sqa; 
                double shu = x(0)*d(0) + x(1)*d(1) + x(2)*d(2);
                double v = 1.0/3.0*shu*sqa;
                area = area + 1.0/2.0*gqcoeff(j)*s; 
                volume = volume + 1.0/2.0*gqcoeff(j)*v;
            }
        }else{ 
            mat ori_dots(A.n_cols-1,3); // ori_dots(11,3);
            for (int j = 1; j < A.n_cols; j++){
                int nodenum = A(i,j);
                ori_dots.row(j-1) = vertex.row(nodenum);
            } 
            mat M(17,11); M.fill(0);  
            mat M1(12,17); M1.fill(0);
            mat M2(12,17); M2.fill(0); 
            mat M3(12,17); M3.fill(0);
            mat M4(11,17); M4.fill(0);
            subdivision_matrix(M, M1, M2, M3, M4);  
            int n = 5;
            //rowvec stotal(3); stotal.fill(0.0); rowvec vtotal(3); vtotal.fill(0.0);  
            for (int j = 0; j < n; j++){
                mat newnodes17 = M*ori_dots; // 17 new nodes
                /////////////////////////////////////////////////
                // element 1
                mat dots = M1*newnodes17; // dots(12,3);
                for (int k =0; k < VWU.n_rows; k++){
                    //rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    mat sf = shape_functions.slice(k);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2); 
                    rowvec d = cross(a_1,a_2)/sqa; 
                    double s = sqa; 
                    double shu = x(0)*d(0) + x(1)*d(1) + x(2)*d(2);
                    double v = 1.0/3.0*shu*sqa;
                    area = area + 1.0/2.0*gqcoeff(k)*s; 
                    volume = volume + 1.0/2.0*gqcoeff(k)*v;
                }
                ///////////////////////////////////////////////
                // element 2 
                dots = M2*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    //rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    mat sf = shape_functions.slice(k);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2); 
                    rowvec d = cross(a_1,a_2)/sqa; 
                    double s = sqa; 
                    double shu = x(0)*d(0) + x(1)*d(1) + x(2)*d(2);
                    double v = 1.0/3.0*shu*sqa;
                    area = area + 1.0/2.0*gqcoeff(k)*s; 
                    volume = volume + 1.0/2.0*gqcoeff(k)*v;
                }
                //////////////////////////////////////////////////////
                // element 3
                dots = M3*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    //rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    mat sf = shape_functions.slice(k);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2); 
                    rowvec d = cross(a_1,a_2)/sqa; 
                    double s = sqa; 
                    double shu = x(0)*d(0) + x(1)*d(1) + x(2)*d(2);
                    double v = 1.0/3.0*shu*sqa;
                    area = area + 1.0/2.0*gqcoeff(k)*s; 
                    volume = volume + 1.0/2.0*gqcoeff(k)*v;
                }
                mat dots4 = M4*newnodes17;    // element 4, still irregular patch
                ori_dots = dots4;
            }
        }
        element_area(i) = area;
        element_volume(i) = volume;
    }
}

double area_constraint_energy_local(mat vertex, Mat<int> face, Mat<int> A, Param param, rowvec elementS_target){
    double E = 0.0;
    double us = param.us;
    int GaussQuadratureN = param.GaussQuadratureN;
    mat VWU = setVMU(GaussQuadratureN); 
    rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN); 
    
    for (int i = 0; i < face.n_rows; i++){
        double s = elementS_target(i);
        if  ( A(i,0) != A(i,1) ){
            mat dots(A.n_cols,3); //dots(12,3); 12 nodes
            for (int j = 0; j < A.n_cols; j++){
                int nodenum = A(i,j);
                dots.row(j) = vertex.row(nodenum);
            }
            // Gaussian quadrature, 3 points
            double a0 = s * 2.0;
            double usi = us;
            for (int j = 0; j < VWU.n_rows; j++){
                rowvec vwu = VWU.row(j); mat sf(12,7); shapefunctions(vwu, sf);
                rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                double  sqa = norm(cross(a_1,a_2),2);  
                E = E + 0.5*usi* 1.0/2.0*gqcoeff(j)*pow(sqa-a0,2.0);
            }
        }else{

            mat ori_dots(A.n_cols-1,3); // ori_dots(11,3);
            for (int j = 1; j < A.n_cols; j++){
                int nodenum = A(i,j);
                ori_dots.row(j-1) = vertex.row(nodenum);
            } 
            mat M(17,11); M.fill(0);  
            mat M1(12,17); M1.fill(0);
            mat M2(12,17); M2.fill(0); 
            mat M3(12,17); M3.fill(0);
            mat M4(11,17); M4.fill(0);
            subdivision_matrix(M, M1, M2, M3, M4);  
            int n = 5;
            //rowvec stotal(3); stotal.fill(0.0); rowvec vtotal(3); vtotal.fill(0.0);  
            for (int j = 0; j < n; j++){
                mat newnodes17 = M*ori_dots; // 17 new nodes
                /////////////////////////////////////////////////
                // element 1
                mat dots = M1*newnodes17; // dots(12,3);
                
                double a0 = s/pow(4.0,j+1.0) * 2.0;
                double usi = us;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                ///////////////////////////////////////////////
                // element 2 
                dots = M2*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                //////////////////////////////////////////////////////
                // element 3
                dots = M3*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                mat dots4 = M4*newnodes17;    // element 4, still irregular patch
                ori_dots = dots4;
            }

        }
    }
    return E;
}

double area_constraint_energy_insertions(mat vertex, Mat<int> face, Mat<int> A, Param param, Mat<int> insertionpatch){
    double E = 0.0;
    double us = param.us;
    int GaussQuadratureN = param.GaussQuadratureN;
    mat VWU = setVMU(GaussQuadratureN); 
    rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN);
    
    for (int m = 0; m < insertionpatch.n_rows; m++){
      for (int n = 0; n < insertionpatch.n_cols; n++){
        int i = insertionpatch(m,n); 
        double s = param.s0;
        if  ( A(i,0) != A(i,1) ){
            mat dots(A.n_cols,3); //dots(12,3); 12 nodes
            for (int j = 0; j < A.n_cols; j++){
                int nodenum = A(i,j);
                dots.row(j) = vertex.row(nodenum);
            }
            // Gaussian quadrature, 3 points
            double a0 = s * 2.0;
            double usi = us;
            for (int j = 0; j < VWU.n_rows; j++){
                rowvec vwu = VWU.row(j); mat sf(12,7); shapefunctions(vwu, sf);
                rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                double  sqa = norm(cross(a_1,a_2),2);  
                E = E + 0.5*usi* 1.0/2.0*gqcoeff(j)*pow(sqa-a0,2.0);
            }
        }else{

            mat ori_dots(A.n_cols-1,3); // ori_dots(11,3);
            for (int j = 1; j < A.n_cols; j++){
                int nodenum = A(i,j);
                ori_dots.row(j-1) = vertex.row(nodenum);
            } 
            mat M(17,11); M.fill(0);  
            mat M1(12,17); M1.fill(0);
            mat M2(12,17); M2.fill(0); 
            mat M3(12,17); M3.fill(0);
            mat M4(11,17); M4.fill(0);
            subdivision_matrix(M, M1, M2, M3, M4);  
            int n = 5;
            //rowvec stotal(3); stotal.fill(0.0); rowvec vtotal(3); vtotal.fill(0.0);  
            for (int j = 0; j < n; j++){
                mat newnodes17 = M*ori_dots; // 17 new nodes
                /////////////////////////////////////////////////
                // element 1
                mat dots = M1*newnodes17; // dots(12,3);
                
                double a0 = s/pow(4.0,j+1.0) * 2.0;
                double usi = us;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                ///////////////////////////////////////////////
                // element 2 
                dots = M2*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                //////////////////////////////////////////////////////
                // element 3
                dots = M3*newnodes17;  
                for (int k =0; k < VWU.n_rows; k++){
                    rowvec vwu = VWU.row(k); mat sf(12,7); shapefunctions(vwu, sf);
                    rowvec x(3); trans_time(sf.col(0),dots,x); //x = strans(sf.col(0)) * dots; //x=sf(:,1)'*dots; 
                    rowvec a_1(3); trans_time(sf.col(1), dots,a_1); //a_1 = strans(sf.col(1)) * dots; //a_1=sf(:,2)'*dots; 
                    rowvec a_2(3); trans_time(sf.col(2),dots,a_2); //a_2 = strans(sf.col(2)) * dots; // a_2=sf(:,3)'*dots;
                    double sqa = norm(cross(a_1,a_2),2);
                    E = E + 0.5*usi* 1.0/2.0*gqcoeff(k)*pow(sqa-a0,2.0);
                }
                mat dots4 = M4*newnodes17;    // element 4, still irregular patch
                ori_dots = dots4;
            }

        }
      }
    }
    return E;
}

vec calculate_all_triangles_area(mat vertex, Mat<int> face){
    vec areas(face.n_rows);
    for (int i = 0; i < face.n_rows; i++){
        int node0 = face(i,0);
        int node1 = face(i,1);
        int node2 = face(i,2);
        rowvec vector0 = vertex.row(node0) - vertex.row(node1); 
        rowvec vector1 = vertex.row(node1) - vertex.row(node2); 
        rowvec vector2 = vertex.row(node2) - vertex.row(node0);
        double side0 = norm(vector0,2.0); 
        double side1 = norm(vector1,2.0); 
        double side2 = norm(vector2,2.0);
        double s = (side0 + side1 + side2)/2.0;
        areas(i) = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
    }
    return areas;
}

rowvec setup_elementS_target(double S0, mat vertex, Mat<int> face, Mat<int> insertionpatch){
    rowvec elementS_target(face.n_rows);
    /*
    // treat the insertions selectedly
    int InsertNumber = insertionpatch.n_rows; // how many insertions
    int PatchNumber = insertionpatch.n_cols;  // one insertion occupies how many patches.
    double insertpatcharea = 2.0 / PatchNumber;     // one insertion's footprint is 2nm2
        double commonpatcharea = (S0 - 2.0*InsertNumber)/( face.n_rows - InsertNumber*PatchNumber );
        for( int i = 0; i < face.n_rows; i ++){
            elementS_target(i) = commonpatcharea;
        }
        for (int i = 0; i < InsertNumber; i++){
            for (int j = 0; j < PatchNumber; j++){
                int facenumber = insertionpatch(i,j);
                elementS_target(facenumber) = insertpatcharea;
            }
        }
        */
    elementS_target.fill(S0/face.n_rows);
    return elementS_target;
} 

double calculate_element_volume(mat dotsin, mat dotsout, int GaussQuadratureN){
    double volume = 0.0;
    mat VWU = setVMU(GaussQuadratureN); 
    rowvec gqcoeff = setVMUcoefficient(GaussQuadratureN);
    for (int i = 0; i < VWU.n_rows; i++) {
        rowvec vwu = VWU.row(i);
        mat sf(12,7);
        shapefunctions(vwu,sf);          // 12 shape functions
        // a_1,2,3 covariant vectors; a1,2 contravariant vectors;
        // inner layer
        rowvec x(3); trans_time(sf.col(0),dotsin,x);
        rowvec a_1(3); trans_time(sf.col(1),dotsin,a_1);
        rowvec a_2(3); trans_time(sf.col(2),dotsin,a_2);
        rowvec xa = cross(a_1,a_2);
        double sqa = norm(xa,2);
        rowvec a_3 = xa/sqa;
        rowvec d = a_3;
        mat shuin = x*strans(d)*sqa;
        // outer layer
        trans_time(sf.col(0),dotsout,x);
        trans_time(sf.col(1),dotsout,a_1);
        trans_time(sf.col(2),dotsout,a_2);
        xa = cross(a_1,a_2);
        sqa = norm(xa,2);
        a_3 = xa/sqa;
        d = a_3;
        mat shuout = x*strans(d)*sqa;
        /////////////////////////////////////////////////////////////
        // volume
        double v = 1.0/3.0*shuout(0,0)-1.0/3.0*shuin(0,0);
        volume = volume + 1.0/2.0*gqcoeff(i)*v;
    }
    return volume;
}

rowvec determine_element_volume(mat vertexin, mat vertexout, Mat<int> face, Mat<int> A, int GaussQuadratureN, double v0, Mat<int> insertionpatch){
    int numface = face.n_rows;
    rowvec elementv(numface);
    // five matrix used for subdivision of the irregular patch
    mat M(17,11); M.fill(0);
    mat M1(12,17); M1.fill(0);
    mat M2(12,17); M2.fill(0);
    mat M3(12,17); M3.fill(0);
    mat M4(11,17); M4.fill(0);
    subdivision_matrix(M, M1, M2, M3, M4);
    for (int i = 0; i < numface; i++){
        double v = 0.0;
        if ( A(i,0) != A(i,1) ) { // for regular patch
            // one ring vertices
            mat dotsin(12,3); 
            mat dotsout(12,3);
            for (int j = 0; j < 12; j++) {
                int nodenum = A(i,j);
                dotsin.row(j) = vertexin.row(nodenum);
                dotsout.row(j) = vertexout.row(nodenum);
            }
            v = calculate_element_volume(dotsin, dotsout, GaussQuadratureN);
        }else{
            int n = 5; // subdivision times
            mat ori_dotsin(11,3);
            mat ori_dotsout(11,3);
            for (int k = 0; k < 11; k++) {
                int nodenum = A(i,k+1);
                ori_dotsin.row(k) = vertexin.row(nodenum);
                ori_dotsout.row(k) = vertexout.row(nodenum);
            }
            mat temp = eye(11,11);
            for (int j = 0; j < n; j++) {
                mat newnodes17in = M*ori_dotsin; // 17 new nodes
                mat newnodes17out = M*ori_dotsout; // 17 new nodes

                if (j != 0) {
                    temp = (M4*M) * temp;
                }
                mat matrix = M*temp;

                mat dotsin = M1*newnodes17in;    // element 1
                mat dotsout = M1*newnodes17out; 
                double v1 = calculate_element_volume(dotsin, dotsout, GaussQuadratureN);
                v = v + v1;

                dotsin = M2*newnodes17in;    // element 2
                dotsout = M2*newnodes17out;
                double v2 = calculate_element_volume(dotsin, dotsout, GaussQuadratureN);
                v = v + v2;

                dotsin = M3*newnodes17in;    // element 3
                dotsout = M3*newnodes17out;
                double v3 = calculate_element_volume(dotsin, dotsout, GaussQuadratureN);
                v = v + v3;

                mat dots4in = M4*newnodes17in;   // element 4, still irregular patch
                mat dots4out = M4*newnodes17out;
                ori_dotsin = dots4in;
                ori_dotsout = dots4out;
            }
        }
        bool isInsertionPatch = false;
        for (int j = 0; j < insertionpatch.n_rows; j++){
            for (int m = 0; m < insertionpatch.n_cols; m++){
                if ( i == insertionpatch(j,m) )
                    isInsertionPatch = true;
            }
        }
        if (isInsertionPatch == true){
            v = v + v0/insertionpatch.n_cols;
        }
        elementv(i) = v;
    }
    return elementv;
}

rowvec calculate_gama(mat vertex, Mat<int> face){
    rowvec gama(face.n_rows);
    for (int i = 0; i < face.n_rows; i++){
            int node0 = face(i,0); // three nodes of this face element
            int node1 = face(i,1);
            int node2 = face(i,2);
            rowvec vector0 = vertex.row(node0) - vertex.row(node1);  double side0 = norm(vector0,2.0); 
            rowvec vector1 = vertex.row(node1) - vertex.row(node2);  double side1 = norm(vector1,2.0); 
            rowvec vector2 = vertex.row(node2) - vertex.row(node0);  double side2 = norm(vector2,2.0);
            double meanside = 1.0/3.0*(side0 + side1+ side2);
            gama(i) = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
    }
    return gama;
}

double calculate_insertionpatch_elementarea(Mat<int> insertionpatch, Mat<double> vertex, Mat<int> face){
    double area = 0.0;
    for (int i = 0; i < insertionpatch.n_rows; i++){
        for (int j = 0; j < insertionpatch.n_cols; j++){
            int facenumber = insertionpatch(i,j);
            int node0 = face(facenumber,0);
            int node1 = face(facenumber,1);
            int node2 = face(facenumber,2);
            rowvec vector0 = vertex.row(node0) - vertex.row(node1); 
            rowvec vector1 = vertex.row(node1) - vertex.row(node2); 
            rowvec vector2 = vertex.row(node2) - vertex.row(node0);
            double side0 = norm(vector0,2.0); 
            double side1 = norm(vector1,2.0); 
            double side2 = norm(vector2,2.0);
            double s = (side0 + side1 + side2)/2.0;
            double thispatcharea = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
            area = area + thispatcharea; 
        }
    }
    double patchnum = insertionpatch.n_rows * insertionpatch.n_cols;
    return area/patchnum;
}

double determine_sideref_global(int i, int j, int facenumber, Mat<int> face, mat vertex, mat vertexold, double gamashu, double meantriS, Mat<int> insertionpatch, bool& isDeformArea, bool& isDeformShape){
    int node0 = face(facenumber,0); // three nodes of this face element
    int node1 = face(facenumber,1);
    int node2 = face(facenumber,2);
    rowvec vector0 = vertex.row(node0) - vertex.row(node1);  double side0 = norm(vector0,2.0); 
    rowvec vector1 = vertex.row(node1) - vertex.row(node2);  double side1 = norm(vector1,2.0); 
    rowvec vector2 = vertex.row(node2) - vertex.row(node0);  double side2 = norm(vector2,2.0);
    double s = (side0 + side1 + side2)/2.0;
    double area = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
    double meanside = 1.0/3.0*(side0 + side1+ side2);
    double gama = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
    if (gama > gamashu){
        isDeformShape = true;
    }
    if ( abs(area-meantriS)/meantriS >= 0.2 ){
        //isDeformArea = true;
    }
    bool isInsertionPatch = false;
    for (int i = 0; i < insertionpatch.n_rows; i++){
        for (int j = 0; j < insertionpatch.n_cols; j++){
            if ( facenumber == insertionpatch(i,j) ){
                isInsertionPatch = true;
            }
        }
    }
    double sideref = 0.0;
    if ( isDeformShape == false && isDeformArea == false ){
        rowvec sideoldvec = vertexold.row(j) - vertexold.row(i); 
        sideref = norm(sideoldvec,2.0);
    }else if ( isDeformArea == true || isInsertionPatch == true){
        sideref = sqrt( 4.0*meantriS / sqrt(3.0) ); 
    }else if ( isDeformShape == true && isDeformArea == false ){
        sideref = sqrt( 4.0*area / sqrt(3.0) ); 
    }
    return sideref;
}

double determine_sideref_local(int i, int j, int facenumber, Mat<int> face, mat vertex, mat vertexold, double gamashu, bool& isDeformShape){
    int node0 = face(facenumber,0); // three nodes of this face element
    int node1 = face(facenumber,1);
    int node2 = face(facenumber,2);
    rowvec vector0 = vertex.row(node0) - vertex.row(node1);  double side0 = norm(vector0,2.0); 
    rowvec vector1 = vertex.row(node1) - vertex.row(node2);  double side1 = norm(vector1,2.0); 
    rowvec vector2 = vertex.row(node2) - vertex.row(node0);  double side2 = norm(vector2,2.0);
    double s = (side0 + side1 + side2)/2.0;
    double area = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
    double meanside = 1.0/3.0*(side0 + side1+ side2);
    double gama = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
    if (gama > gamashu){
        isDeformShape = true;
    }
    double sideref = 0.0;
    if ( isDeformShape == false ){
        rowvec sideoldvec = vertexold.row(j) - vertexold.row(i); 
        sideref = norm(sideoldvec,2.0);
    }else if ( isDeformShape == true ){
        sideref = sqrt( 4.0*area / sqrt(3.0) ); 
    }
    return sideref;
}

rowvec determine_spontaneous_curvature(Param param, Mat<int> face, mat vertex, Mat<int> insertionpatch){
    double C0 = param.C0;
    double c0 = param.c0;
    double R  = param.R;
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
            double dismin = 2.0*R;
            for (int k = 0; k < insertionpatch.n_cols; k++){
                int facenumber = insertionpatch(j,k);
                rowvec center0 = 1.0/3.0 * ( vertex.row(face(facenumber,0)) + vertex.row(face(facenumber,1)) + vertex.row(face(facenumber,2)) );
                rowvec disvec = centeri - center0;
                double dis = norm(disvec,2.0);
                double distance = 2.0*R * asin(dis/2.0/R);
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
            double dismin = 2.0*R;
            for (int k = 0; k < insertionpatch.n_cols; k++){
                int facenumber = insertionpatch(j,k);
                rowvec center0 = 1.0/3.0 * ( vertex.row(face(facenumber,0)) + vertex.row(face(facenumber,1)) + vertex.row(face(facenumber,2)) );
                rowvec disvec = centeri - center0;
                double dis = norm(disvec,2.0);
                double distance = 2.0*R * asin(dis/2.0/R);
                if (distance < dismin){
                    dismin = distance;
                }
            }
            spontcurvs(j) = abs(C0) * exp( - pow(dismin/sigma,2.0)/2.0 );
        }
        if (isAdditiveScheme == true){
            double shu = 0.0;
            double H = 1.0/R;
            for (int j = 0; j < insertionpatch.n_rows; j++){
                shu = shu + pow(2.0*H-spontcurvs(j),2.0);
            }
            shu = shu - pow(2.0*H,2.0);
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
        }else{
            spontcurv(i) = - max(spontcurvs); // here the spont_curvature is negative sign
        }
        
        if ( abs(spontcurv(i)) < 1.0e-15 )
            spontcurv(i) = 0.0;
    }
    return spontcurv;
}

mat update_vertex(mat vertex, double a, mat force, Mat<int> insertionpatch, Mat<int> face,Param param){
    mat vertexnew = vertex + a*force;
    return vertexnew;
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
