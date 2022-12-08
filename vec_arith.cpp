#include "CJVector.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------VECTOR carithmetic functions-----------------------------------------------//
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

//Elementwise Addition
CJVector operator+(const CJVector& lhs, const CJVector& rhs){
    if(!dimMatch(lhs,rhs))
        return CJVector();
    CJVector res(lhs.len);
    for(size_t i=0; i<lhs.len; i++)
        res[i]=lhs[i]+rhs[i];
    return res;
}

//Elementwise Subtraction
CJVector operator-(const CJVector& lhs, const CJVector& rhs){
    if(!dimMatch(lhs,rhs))
        return CJVector();
    CJVector res(lhs.len);
    for(size_t i=0; i<lhs.len; i++)
        res[i]=lhs[i]-rhs[i];
    return res;
}