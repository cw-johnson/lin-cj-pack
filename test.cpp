#include <iostream> //I/O
#include <chrono>
#include <random>
#include "CJVector.h"

#if defined(__linux__)
    #include <cblas.h> //Generic BLAS Routines 
#elif __APPLE__
    #include <Accelerate/Accelerate.h> //Apple Optimized BLAS routines
#else   
    # error "Unknown OS"
#endif

using namespace std;


int main(){
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> trad_time, cblas_time;

    uniform_real_distribution<double> myrand(-10.0, 10.0);
    default_random_engine re;
    

    const size_t n = 10;
    double arr1[n] = {1,2,3,4,5,6,7,8,9,10};
    double arr2[n] = {1,2,3,4,5,6,7,8,9,10};
    double res_t = 0;
    double res_o = 0;
    
    start = chrono::system_clock::now();
    end = chrono::system_clock::now();
    trad_time = end - start;

    CJVector A;

    cout<<"Done"<<endl;
    return 0;




}