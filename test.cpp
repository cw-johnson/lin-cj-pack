#include <iostream> //I/O
#include <chrono>
#include <random>
#include "CJVector.h"

#if defined(__linux__)
    #include <cblas.h> //Generic BLAS Routines 
#elif __APPLE__
    #include <Accelerate/Accelerate.h> //Apple Optimized BLAS routines
#else   
    # error "Unknown Compiler"
#endif

using namespace std;


int main(){
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> trad_time, cblas_time;

    uniform_real_distribution<double> myrand(-10.0, 10.0);
    default_random_engine re;
    


    start = chrono::system_clock::now();
    end = chrono::system_clock::now();
    trad_time = end - start;

    CJVector A;

    cout<<"Done"<<endl;
    return 0;




}