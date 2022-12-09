#include "CJVector.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------VECTOR arithmetic functions-----------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Math Functions
double dot(const CJVector& lhs, const CJVector& rhs){
    if (!dimMatch(lhs,rhs))
        return NAN;
    return cblas_ddot(lhs.len,lhs.data,1,rhs.data,1);
}

//Euclidean Norm
double norm(const CJVector& arg){
    return cblas_dnrm2(arg.len,arg.data,1);
}

//Vector Addition
CJVector operator+(const CJVector& lhs, const CJVector& rhs){
    if(!dimMatch(lhs,rhs))
        return CJVector();
    CJVector res(lhs.len);
    for(size_t i=0; i<lhs.len; i++)
        res[i]=lhs[i]+rhs[i];
    return res;
}

//Vector Subtraction
CJVector operator-(const CJVector& lhs, const CJVector& rhs){
    if(!dimMatch(lhs,rhs))
        return CJVector();
    CJVector res(lhs.len);
    for(size_t i=0; i<lhs.len; i++)
        res[i]=lhs[i]-rhs[i];
    return res;
}

//Scalar-Vector Multiplication
CJVector operator*(const double &a, const CJVector& x){
    CJVector res = x;
    cblas_dscal(x.len, a, res.data, 1);
    return res;
}