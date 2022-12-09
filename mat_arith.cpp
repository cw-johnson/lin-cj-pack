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
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.m, rhs.n, lhs.n, 1, lhs.data, lhs.m, rhs.data, rhs.m, 0, res.data, 1);
    return res;
}