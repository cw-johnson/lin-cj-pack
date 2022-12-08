/*
 * Library:      CJ's Linear Algebra Library V2
 * Course:       MAD 5420: Numerical Optimization, Fall '22
 * Instructor:   Professor Dr. Kyle Gallivan
 * Student Name: Christopher Johnson
 * Date:         12/07/2022
 * Description:  
*/



#include <iostream> //I/O
#include<iomanip> //Output formatting
#include <chrono> //Measure Timing
#include <cstring> //memcpy instruction
#include <cmath> //Math Stuff
#include <initializer_list> //Initializer List Construction

//Detect OS
#if defined(__linux__)
    #include <cblas.h> //Generic BLAS Routines 
#elif __APPLE__
    #include <Accelerate/Accelerate.h> //Apple Optimized BLAS routines
#else   
    # error "Unknown OS"
#endif

using namespace std;

#ifndef CJVECTOR_H
#define CJVECTOR_H

class CJVector;
class CJMatrix; //Forward Declarations

class CJVector{
    public:
        //Default Constructor
        CJVector() : len{0}, trans{false}, data{nullptr} {} 
        //Parameterized Constructor
        CJVector(size_t v_len) : len{v_len}, trans{false}, data{new double[v_len]} {
            zeros();
        } 
        //Initializer List Constructor
        CJVector(initializer_list<double> list) : len(list.size()), trans{false}, data{new double[len]} {
            size_t count = 0;
            for (auto i : list)
                data[count++] = i;
        }
        //Destructor
        ~CJVector(){
            if (data != nullptr)
                delete [] data; 
            data=nullptr; len=0;
        }
        //Copy Constructor
        CJVector(const CJVector &rhs){ 
            this->len = rhs.len; 
            this->trans = rhs.trans;
            data = new double[len];
            memcpy(data, rhs.data, len * 8);
        }
        //Move Constructor
        CJVector(CJVector&& rhs){ 
            this->len = rhs.len;
            this->trans = rhs.trans;
            this->data = rhs.data;
            rhs.data = nullptr;
            rhs.len = 0;
        }
        //Copy Assignment Operator
        CJVector& operator=(const CJVector& rhs){
            if(this != &rhs){
                if (this->len != rhs.len || sizeof(data) != sizeof(rhs.data)){
                    if (this->data != nullptr)
                        delete[] this->data;
                    if (rhs.len != 0)
                        this->data = new double[rhs.len];
                }
                this->len = rhs.len;
                this->trans = rhs.trans;
                memcpy(this->data, rhs.data, this->len * 8);
            }
            return *this;
        }


        //Vector Arithmetic Functions
        friend double dot(const CJVector& lhs, const CJVector& rhs); //Dot Product O(n)
        friend double norm(const CJVector& arg); //Euclidean Norm O(n)
        friend CJVector operator+(const CJVector& lhs, const CJVector& rhs); //Elementwise Addition O(n)
        friend CJVector operator-(const CJVector& lhs, const CJVector& rhs); //Elementwise Subtraction O(n)
        friend CJVector operator*(const double &a, const CJVector& x); //Scalar Multiplication

        //Shared Arithmetic Functions
        friend CJVector operator*(const CJMatrix &lhs, const CJVector &rhs);

        //Utility
        double& operator[] (const size_t& i){return data[i];}
        const double& operator[](const size_t& i) const {return data[i];}
        void print();
        void setAll(double val);
        void zeros();
        void addRandomData(double min, double max);
        friend bool dimMatch(const CJVector& lhs, const CJVector& rhs);
    private:
        size_t len;
        bool trans;
        double * data; 
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////






enum MatrixType{DENSE,SYMMETRIC};

class CJMatrix{
    public:
        //Default Constructor
        CJMatrix() : len{0}, m{0}, n{0}, trans{false}, mtype{DENSE}, data{nullptr} {} 
        //Parameterized Constructor
        CJMatrix(size_t m, size_t n) : len{m*n}, m{m}, n{n}, trans{false}, mtype{DENSE}, data{new double[m*n]} {
            zeros();
        } 
        CJMatrix(size_t m, size_t n, MatrixType tt) : len{m*n}, m{m}, n{n}, trans{false}, mtype{tt}, data{new double[m*n]} {
            zeros();
        }         
        //Initializer List Constructor
        // CJMatrix(initializer_list<double> list) : len(list.size()), trans{false},  mtype{DENSE}, data{new double[len]} {
        //     size_t count = 0;
        //     for (auto i : list)
        //         data[count++] = i;
        // }
        //Destructor
        ~CJMatrix(){
            if (data != nullptr)
                delete [] data; 
            data=nullptr; len=0; m=0; n=0; trans=false;
        }
        //Copy Constructor
        CJMatrix(const CJMatrix &rhs){ 
            this->len = rhs.len; 
            this->m = rhs.m;
            this->n = rhs.n;
            this->trans = rhs.trans;
            this->mtype = rhs.mtype;
            data = new double[len];
            memcpy(data, rhs.data, len * 8);
        }
        //Move Constructor
        CJMatrix(CJMatrix&& rhs){ 
            this->len = rhs.len;
            this->m = rhs.m;
            this->n = rhs.n;
            this->trans = rhs.trans;
            this->mtype = rhs.mtype;
            this->data = rhs.data;
            rhs.data = nullptr;
            rhs.len = 0;
        }
        //Copy Assignment Operator
        CJMatrix& operator=(const CJMatrix& rhs){
            if(this != &rhs){
                if (this->len != rhs.len || sizeof(data) != sizeof(rhs.data)){
                    if (this->data != nullptr)
                        delete[] this->data;
                    if (rhs.len != 0)
                        this->data = new double[rhs.len];
                }
                this->len = rhs.len;
                this->m = rhs.m;
                this->n = rhs.n;
                this->mtype = rhs.mtype;
                this->trans = rhs.trans;
                memcpy(this->data, rhs.data, this->len * 8);
            }
            return *this;
        }


        //Row Access
        double& operator[] (const int& row);
        const double &operator[](const int &row) const;
        //Element Access
        double& operator() (const int& row, const int& col);
        const double& operator() (const int& row,const int& col) const;
        double& getelem(const int& row, const int& col);
        const double& getelem(const int& row,const int& col) const;
        //Utility
        int calcIndex(const int row, const int col) const;
        void print();
        friend bool dimMatch(const CJMatrix& lhs, const CJMatrix& rhs);
        //Matrix Arithmetic


        //Shared Arithmetic
        friend CJVector operator*(const CJMatrix &lhs, const CJVector &rhs);


        void setAll(double val);
        void zeros();
        void addRandomData(double min, double max);
    private:
        size_t len, m, n; //
        bool trans;
        MatrixType mtype;
        double * data; 
};



#include "mat_core.cpp" //Matrix Utility Functions
#include "vec_core.cpp" //Vector Utility Functions

#include "mat_arith.cpp" //Matrix Arithmetic Functions
#include "vec_arith.cpp" //Vector Arithmetic Functions



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------COMBINED MATRIX/VECTOR FUNCTIONS--------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Optimized Matrix-Vector Product (Support for Symmetric Matricies)
CJVector operator*(const CJMatrix &lhs, const CJVector &rhs){
    CJVector res(lhs.m);
    if(lhs.mtype==SYMMETRIC)
        cblas_dsymv(CblasRowMajor, CblasUpper, lhs.n, 1, lhs.data, lhs.n, rhs.data, 1, 0, res.data, 1);
    else
        cblas_dgemv(CblasRowMajor, CblasNoTrans, lhs.m, lhs.n, 1, lhs.data, lhs.m, rhs.data, 1, 0, res.data, 1);
    return res;

}








#endif
