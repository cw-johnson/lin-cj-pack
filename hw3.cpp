
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



struct OptResult{
    string method = "N/A"; //Optimization Method
    string problem_type = "N/A"; //Problem Type
    unsigned int dimension = 0; //Problem in R^n
    string problem_ID = "N/A"; //Problem ID
    unsigned int num_iterations = 0; //Number of iterations until convergence
    unsigned int num_l1_blas = 0; //# of level 1 BLAS routines used
    unsigned int num_l2_blas = 0; //# of level 2 BLAS routines used
    unsigned int num_l3_blas = 0; //# of level 3 BLAS routines used

    chrono::duration<double> ex_time = chrono::duration<double>::zero(); //Measured Execution Time
    friend ostream& operator<<(ostream& os, const OptResult &res){
        os<< setw(PRINT_W) << left << "Method:"              << res.method               <<endl;
        os<< setw(PRINT_W) << left << "Problem Type:"        << res.problem_type         <<endl;
        os<< setw(PRINT_W) << left << "Problem Dimension:"   << res.dimension            <<endl;
        os<< setw(PRINT_W) << left << "Problem ID:"          << res.problem_ID           <<endl;
        os<< setw(PRINT_W) << left << "# Iterations:"        << res.num_iterations       <<endl;
        os<< setw(PRINT_W) << left << "Execution Time (ms):" << 1000*res.ex_time.count() <<endl;
        os<< setw(PRINT_W) << left << "# Level 1 BLAS:"      << res.num_l1_blas          <<endl;
        os<< setw(PRINT_W) << left << "# Level 2 BLAS:"      << res.num_l2_blas          <<endl;
        os<< setw(PRINT_W) << left << "# Level 3 BLAS:"      << res.num_l3_blas          <<endl;
        return os;
    }
};


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
    test.num_l1_blas = 1;
    test.num_l2_blas = 2;
    test.num_l3_blas = 3;
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

    CJVector A = {2.0, 2.0};
    A.print(); 

    CJVector B = {3.0, 3.0}; 
    B.print();

    double d = dot(A,B);
    cout<<d<<endl;

    d = norm(A+B);
    cout<<d<<endl;

    CJVector D = A;
    
    testBLAS();
    return 0;
}