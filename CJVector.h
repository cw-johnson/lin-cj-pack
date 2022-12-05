#include <iostream> //I/O
#include <chrono> //Measure Timing
#include <cstring> //memcpy instruction

#if defined(__linux__)
    #include <cblas.h> //Generic BLAS Routines 
#elif __APPLE__
    #include <Accelerate/Accelerate.h> //Apple Optimized BLAS routines
#else   
    # error "Unknown Compiler"
#endif

using namespace std;

#ifndef CJVECTOR_H
#define CJVECTOR_H



class CJVector{
    public:
        //Default Constructor
        CJVector() : len{0}, trans{false}, data{nullptr} {} 
        //Parameterized Constructor
        CJVector(unsigned int v_len) : len{v_len}, trans{false}, data{new double[v_len]} {} 
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
        friend CJVector operator+(const CJVector& lhs, const CJVector& rhs);
        friend CJVector operator-(const CJVector& lhs, const CJVector& rhs);
        friend CJVector dot(const CJVector& lhs, const CJVector& rhs);
    private:
        size_t len;
        bool trans;
        double * data; 
};

#endif

CJVector operator+(const CJVector& lhs, const CJVector& rhs){
    CJVector res(lhs.len);
    for (size_t i = 0; i<lhs.len; i++)
        res.data[i] = lhs.data[i];
    return res; 
}