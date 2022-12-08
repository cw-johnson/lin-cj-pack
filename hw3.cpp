
/*
 * Assignment:   Graded HW #3: Inexact Newton, Quasi-Newton, & Steepest Descent
 * Course:       MAD 5420: Numerical Optimization, Fall '22
 * Instructor:   Professor Dr. Kyle Gallivan
 * Student Name: Christopher Johnson
 * Date:         12/07/2022
 * Description:  
*/


#include <iostream> //...for Basic I/O
#include <iomanip> //.. for Output formatting
#include <chrono> //.. for Execution Time measurement
#include <cmath> //.. for pow()
#include <random> //.. for producing random values
#include <limits> //.. for infinity(), unit roundoff
#include <algorithm> //.. for max()
#include <string.h> //.. for Naming things
#include <vector> //Memory management
#include "CJVector.h" //Legacy Matrix Implementation

//Detect OS
#if defined(__linux__)
    #include <cblas.h> //Generic BLAS Routines 
#elif __APPLE__
    #include <Accelerate/Accelerate.h> //Apple Optimized BLAS routines
#else   
    # error "Unknown OS"
#endif


using namespace std;

#define PRINT_W 25


/*
 * CBLAS notes
 * D prefix -> for Doubles
 * S -> Singles
 * C -> Complex
 * Z -> Complex16
 * -> We only need D prefix functions
 * 
 * L1: O(n)
 * Dot Product: cblas_ddot
 * Scale:       cblas_dscal
 * MAC:         cblas_daxpy
 * L2 Norm      cblas_dnrm2
 * Copy         cblas_dcopy
 * Swap         cblas_dswap
 * 
 * L2: O(n^2)
 * Matrix-Vector Product:           cblas_dgemv
 * Symmetric Matrix-Vector Product: cblas_dsymv
 * 
 * L3: O(n^3)
 * Matrix-Matrix Product:           cblas_dgemm
 * Symmetrix Matrix-Matrix Product: cblas_dsymm
*/


//Unconstrained Quadratic Programming Problem 
struct QP_Problem{
    string name;
    size_t n;
    CJMatrix A; //A must be s.p.d
    CJVector b;
    double f(const CJVector &x){
        double res = dot(x,A*x); //n^2 + n
        res += dot(b,x); //2n
        return res;
    }
    CJVector grad(const CJVector &x){
        return A*x+b; //n^2 + n
    }
    CJVector r(const CJVector &x){
        return b-A*x; //n^2 + n
    }
};


//Steepest Descent Method Update
//Operations: 2n^2 + 4n
unsigned long int SD_Update(QP_Problem &qp, CJVector &xk){
    CJVector pk = qp.r(xk); //n^2 + n
    double ak = dot(pk,pk)/dot(pk,qp.A*pk); //n^2 + 2n 
    xk = xk + ak*pk; //2n
    unsigned long int n = qp.n;
    return 2*pow(n,2) + 4*n;
}

//Richardson's Stationary Method Update
//Operations: n^2 + 3n
unsigned long int RS_Update(QP_Problem &qp, CJVector &xk, const double &ak){
    CJVector pk = qp.r(xk); //n^2 + n;
    xk = xk + ak*pk; //2n
    unsigned long int n = qp.n;
    return pow(n,2) + 3*n;
}



struct OptResult{
    string method = "N/A"; //Optimization Method
    string problem_type = "N/A"; //Problem Type
    unsigned int dimension = 0; //Problem in R^n
    string problem_ID = "N/A"; //Problem ID
    unsigned int num_iterations = 0; //Number of iterations until convergence
    unsigned long int num_ops = 0; //# of elementary operations used


    chrono::duration<double> ex_time = chrono::duration<double>::zero(); //Measured Execution Time
    friend ostream& operator<<(ostream& os, const OptResult &res){
        os<< setw(PRINT_W) << left << "Method:"              << res.method              <<endl;
        os<< setw(PRINT_W) << left << "Problem Type:"        << res.problem_type        <<endl;
        os<< setw(PRINT_W) << left << "Problem Dimension:"   << res.dimension           <<endl;
        os<< setw(PRINT_W) << left << "Problem ID:"          << res.problem_ID          <<endl;
        os<< setw(PRINT_W) << left << "# Iterations:"        << res.num_iterations      <<endl;
        os<< setw(PRINT_W) << left << "Execution Time (ms):" << 1000*res.ex_time.count()<<endl;
        os<< setw(PRINT_W) << left << "# Operations:"        << res.num_ops             <<endl;
        return os;
    }
};


OptResult solve_QP(QP_Problem qp, CJVector &xko, string method){
    unsigned int k=0; 
    unsigned int kmax=

}



void testResultStruct(){
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> ex_time;
    start = chrono::system_clock::now();
    end = chrono::system_clock::now();
    ex_time = end-start;
    OptResult test;
    test.method = "Steepest Descent";
    test.problem_type = "Quadratic";
    test.dimension = 3;
    test.problem_ID = "#1, Convex";
    test.num_iterations = 100;
    test.num_ops = 1;
    test.ex_time = ex_time;
    cout<<test<<endl;
}

void testBLAS(){
    const size_t n = 10; const unsigned int inc = 1;
    
    double arr1[n] = {1,2,3,4,5,6,7,8,9,10};
    double arr2[n] = {1,2,3,4,5,6,7,8,9,10};
    double res_t = 0;
    double res_o = 0;
    for(size_t i = 0; i < n; i++){
        res_t += arr1[i]*arr2[i];
    }

    res_o = cblas_ddot(n,arr1,inc,arr2,inc);
    cout<<"Res_T="<<res_t<<endl;
    cout<<"Res_O="<<res_o<<endl;   
}

int main(){
    //Output Header
    cout<<endl;
    cout<<"Graded HW #3: Inexact Newton, Quasi-Newton, & Steepest Descent"<<endl;
    cout<<"MAD 5420: Numerical Optimization, Fall '22"<<endl;
    cout<<"Professor Dr. Kyle Gallivan"<<endl;
    cout<<"Implementation by Christopher Johnson"<<endl;
    cout<<"December 7, 2022"<<endl;
    cout<<endl;

    CJMatrix G(2,2,SYMMETRIC);
    G(0,0)=16.0; G(0,1) = 0.0;
    G(1,0)=0.0;  G(1,1) = 4.0;

    CJVector c = {2.0, 2.0};

    CJVector x = {1.5, 1.5};
    QP_Problem q1 = {"Hello",2,G,c};
    cout<<"q1.f()"<<q1.f(x)<<endl;
    CJVector grad = q1.grad(x);
    grad.print();
    return 0;
}