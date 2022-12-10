#include "CJVector.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------Matrix arithmetic functions-----------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Scalar-Matrix Multiplication
CJMatrix operator*(const double &a, const CJMatrix& x){
    CJMatrix res = x;
    cblas_dscal(x.len, a, res.data, 1);
    return res;
}

//Matrix Addition
CJMatrix operator+(const CJMatrix& lhs, const CJMatrix& rhs){
    if(!dimMatch(lhs,rhs))
        return CJMatrix();
    CJMatrix res(lhs.m, lhs.n);
    for(size_t i=0; i<lhs.m; i++)
        for(size_t j=0; j<lhs.n; j++)
            res(i,j) = lhs(i,j) + rhs(i,j);
    return res;
}

//Matrix Subtraction
CJMatrix operator-(const CJMatrix& lhs, const CJMatrix& rhs){
    if(!dimMatch(lhs,rhs))
        return CJMatrix();
    CJMatrix res(lhs.m,lhs.n);
    for(size_t i=0; i<lhs.m; i++)
        for(size_t j=0; j<lhs.n; j++)
            res(i,j) = lhs(i,j) - rhs(i,j);
    return res;
}

//Matrix Multiplication
CJMatrix operator*(const CJMatrix& lhs, const CJMatrix& rhs){
    if(lhs.n != rhs.m)
        return CJMatrix();
    CJMatrix res(lhs.m,rhs.n);
    if (lhs.mtype == SYMMETRIC)
        cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, lhs.m, rhs.n , 1, lhs.data, lhs.m, rhs.data, lhs.m, 0, res.data, lhs.m); //Symmetric Matrix Support (n^2)
    else
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.m, rhs.n, lhs.n, 1, lhs.data, lhs.m, rhs.data, lhs.n, 0, res.data, lhs.m); //Dense Matrix (n^3)
    return res;
}

//2x2 Symetric Matrix Norm
//Closed form from Deledalle, Denis, Taenable_backtracki, & Tupin
double m2norm_2(CJMatrix &G){
    if((G.m != 2) || (G.n != 2))
        return NAN;
    double a = G(0,0);
    double b = G(1,1);
    double c = G(1,0);
    double delta = sqrt(4*pow(c,2) + pow(a-b,2));
    double l1 = (a+b-delta)/2.0;
    double l2 = (a+b+delta)/2.0;
    return max(abs(l1),abs(l2));
}

//3x3 Symetric Matrix Norm
//Closed form from Deledalle, Denis, Taenable_backtracki, & Tupin
double m2norm_3(CJMatrix &G){
    if((G.m != 3) || (G.n != 3))
        return NAN;
    double a = G(0,0);
    double b = G(1,1);
    double c = G(2,2);
    double d = G(1,0);
    double e = G(2,1);
    double f = G(2,0);
    double x1 = a*a + b*b + c*c - a*b - a*c - b*c + 3*(d*d + f*f + e*e);
    double x2 = -(2*a-b-c)*(2*b-a-c)*(2*c-a-b) +9*((2*c-a-b)*d*d + (2*b-a-c)*f*f + (2*a-b-c)*e*e) - 54*(d*e*f);
    double phi;
    if(x2>0)
        phi = atan(sqrt(4*x1*x1*x1 - x2*x2)/x2);
    else if(x2<0)
        phi = atan(sqrt(4*x1*x1*x1 - x2*x2)/x2) + M_PI;
    else
        phi = M_PI/2;
    double l1 = (a+b+c - 2*sqrt(x1)*cos(phi/3))/3;
    double l2 = (a+b+c + 2*sqrt(x1)*cos((phi-M_PI)/3))/3;
    double l3 = (a+b+c + 2*sqrt(x1)*cos((phi+M_PI)/3))/3;
    return max({abs(l1),abs(l2),abs(l3)});
}