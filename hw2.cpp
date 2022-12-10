
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
#include <string.h> //.. for naming stuff
#include "CJVector.h" //Matrix & Vector Implementation V2

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

//Unconstrained Quadratic Programming Problem 
struct QP_Problem{
    string name;
    size_t n;
    CJMatrix A; //Gradient and H for s.p.d A
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
    CJMatrix H(){
        return 2*A;
    }
};


//Steepest Descent Method Update
//Operations: 2n^2 + 4n
unsigned long int SD_Update(QP_Problem &qp, CJVector &xk){
    CJVector pk = qp.r(xk); //n^2 + n
    double ak = dot(pk,pk)/dot(pk,qp.A*pk); //n^2 + 2n 
    //Take Step
    xk = xk + ak*pk; //2n
    unsigned long int n = qp.n;
    return 2*pow(n,2) + 4*n;
}

//Richardson's Stationary Method Update
//Operations: n^2 + 3n
unsigned long int RS_Update(QP_Problem &qp, CJVector &xk, const double &ak){
    CJVector pk = qp.r(xk); //n^2 + n;
    //Take Step
    xk = xk + ak*pk; //2n
    unsigned long int n = qp.n;
    return pow(n,2) + 3*n;
}

//SDSlow Method Update
//Operations: n^2 + 3n + 1
unsigned long int SDS_Update(QP_Problem &qp, CJVector &xk, double &ak, const double &sigma){
    CJVector pk = qp.r(xk); //n^2 + n;
    ak = sigma*ak; //1
    //Take Step
    xk = xk + ak*pk; //2n
    unsigned long int n = qp.n;
    return pow(n,2) + 3*n + 1;
}

//Conjugate Gradient Method Update
//Operations: n^2 + 8n + 2
unsigned long int CG_Update(QP_Problem &qp, CJVector &xk, CJVector &rk, CJVector &dk, double &sigk){
    CJVector vk = qp.A*dk; //n^2
    double uk = dot(dk,vk); //n
    double ak = sigk/uk; //1
    //Take Step
    xk = xk + ak*dk; //2n
    rk = rk - ak*vk; //2n
    double sig_old = sigk;
    sigk = dot(rk,rk); //n
    double Bk = sigk/sig_old; //1
    dk = rk + Bk*dk; //2n
    unsigned long int n = qp.n;
    return n^2 + 8*n + 2;
}

//Gauss-Southwell Method Update
//Operations: n^2 + 5n operations
unsigned long int GS_Update(QP_Problem &qp, CJVector &xk, const double &ak){
    CJVector rk = qp.r(xk); //n^2 + n;
    double max = 0;
    int maxi = 0;
    for (int i =0; i<qp.n; i++){ //n
        if (abs(rk[i])>max){
        max = abs(rk[i]);
        maxi = i;
        }
    }
    CJVector ei = qp.A.slice(0,maxi,qp.n);
    double dp = dot(ei,rk); //n

    if (dp<0.0) //n
        rk = -1.0*ei; 
    else
        rk = ei;

    //Take Step
    xk = xk + ak*rk; //n

    unsigned long int n = qp.n;
    return pow(n,2) + 5*n;
}




// double backtracking_armijo(double a, double L){
//     double beta = 0.95;
//     double alpha = a; double alpha_new = alpha;
//     double lambda = 0.01; 
//     int m = 0; int mmax = 100;
//     while(m<mmax){
//         if(alpha_new <= (2*(1-lambda))/L){
//             return alpha_new; break;
//         }
//         alpha_new = beta*alpha_new;
//         m++;
//     }
//     //cout<<m<<' '<<alpha_new<<endl;
//     return alpha_new;
// }

double backtracking_wolfe(QP_Problem &qp, CJVector &xk, CJVector &pk, double a0, double beta){
    double a = a0;
    int m = 0; int mmax = 10;
    double lhs1 = 1.0;
    double rhs1 = 0.0;
    double lhs2 = 0.0;
    double rhs2 = 1.0;
    double c1 = 0.25; double c2 = 0.75;
    CJVector gradient = qp.grad(xk);
    //cout<<"Gradient: ";
    //gradient.print();
    while((m<mmax) && ((lhs1>rhs1) || (lhs2<rhs2))){
        lhs1 = qp.f(xk + a*pk); //Sufficient Decrease
        rhs1 = qp.f(xk) + (c1*a)*dot(gradient,pk); 

        // lhs2 = dot(qp.grad(xk + a*pk),pk);
        // rhs2 = c2*dot(qp.grad(xk),pk);

        lhs2 = abs(dot(qp.grad(xk + a*pk),pk));
        rhs2 = (c2*dot(qp.grad(xk),pk));
        //cout << lhs << ' ' << rhs <<endl;
        a = beta*a;
        m++;
    }
    //cout<<m<<' '<<a<<endl;
    return a;
}
unsigned long int NTCG_Update(QP_Problem &qp, CJVector &xk, const double &ak){

    return 0;
}


unsigned long int BFGS_Update(QP_Problem &qp, CJVector &xk, CJMatrix &Hk){
    Hk.print(); cout<<endl;
    CJVector pk = -1.0 * Hk *qp.grad(xk);
    pk.print(); cout<<endl;
    double a0_bt = 1.0; //Backtrack with a0 = 1
    double alpha = backtracking_wolfe(qp,xk,pk,a0_bt,0.99);
    CJVector xk_old = xk;
    xk = xk + alpha*pk;
    cout<<"x: "; xk.print(); 
    cout<<"f(x): " << qp.f(xk) <<endl;
    CJVector sk = xk - xk_old;
    CJVector yk = qp.grad(xk) - qp.grad(xk_old);
    double rhok = 1/dot(yk,sk);
    //CJVector rk = sk - Hk*yk;
    //double muk = dot(rk,yk)/( pow(rhok,2));
    //Hk = Hk + rhok*(outer(rk,sk)+outer(sk,rk)) - muk*outer(sk,sk);
    CJMatrix I(qp.n,qp.n);
    I.makeIdentity();

    Hk = (I-rhok*outer(sk,yk))*Hk*(I-rhok*outer(yk,sk)) + rhok*outer(sk,sk);

    //Hk.print();
    //cout<<endl; cout<<endl;
    return 0;
}
struct OptResult{
    string method = "N/A"; //Optimization Method
    string problem_type = "N/A"; //Problem Type
    unsigned int dimension = 0; //Problem in R^n
    string problem_ID = "N/A"; //Problem ID
    unsigned int num_iterations = 0; //Number of iterations until convergence
    unsigned long int num_ops = 0; //# of elementary operations used
    CJVector solution = CJVector();
    double sigma = 0;
    double lambda_max = 0;
    double f_eval = 0;
    chrono::duration<double> ex_time = chrono::duration<double>::zero(); //Measured Execution Time
    friend ostream& operator<<(ostream& os, const OptResult &res){
        os<< setw(PRINT_W) << left << "Method:"              << res.method              <<endl;
        if(res.method == "SDslow")
            os<< setw(PRINT_W) << left << "Sigma:"           << res.sigma               <<endl;
        else if(res.method == "Richardson's Stationary")
            os<< setw(PRINT_W) << left << "Lambda_Max:"      << res.lambda_max          <<endl;
        os<< setw(PRINT_W) << left << "Problem Type:"        << res.problem_type        <<endl;
        os<< setw(PRINT_W) << left << "Problem Dimension:"   << res.dimension           <<endl;
        os<< setw(PRINT_W) << left << "Problem ID:"          << res.problem_ID          <<endl;
        os<< setw(PRINT_W) << left << "# Iterations:"        << res.num_iterations      <<endl;
        os<< setw(PRINT_W) << left << "Execution Time (ms):" << 1000*res.ex_time.count()<<endl;
        os<< setw(PRINT_W) << left << "# Operations:"        << res.num_ops             <<endl;
        os<< setw(PRINT_W) << left << "f(x*) =      "        << res.f_eval              <<endl;
        os<< setw(PRINT_W) << left << "Solution x* ="        << endl                    <<endl;
        os<<res.solution<<endl;
        return os;
    }
};

struct CGResult{
    CJVector x = CJVector();
    unsigned int num_ops = 0;
};

CGResult solve_modCG(QP_Problem qp, CJVector &xko, double eps){
    unsigned int k=0;
    unsigned int kmax = 100;
    unsigned int op_ctr = 0; 
    CJVector xk = xko;
    CJVector r_cg = qp.r(xk); //Initial Calculations for CG
    CJVector d_cg = r_cg;  //for CG
    double sig_cg = dot(r_cg,r_cg); //..for CG
    while((k<kmax) && (norm(qp.r(xk))>eps)){
        op_ctr += CG_Update(qp,xk,r_cg,d_cg,sig_cg);
        k++;
    }
    CGResult res;
    res.x = xk;
    res.num_ops = op_ctr;
    return res;
}


OptResult solve_QP(QP_Problem &qp, CJVector &xko, string method, double eps){
    chrono::time_point<chrono::system_clock> start, end; //Setup Timing
    chrono::duration<double> ex_time;
    unsigned int op_ctr = 0;    

    unsigned int k=0; 
    unsigned int kmax=10000; //Max no Iterations
    CJVector xk = xko;

    double lambda_max = 0.1; //For Richardsons Stationary
    double a_rs = 1/lambda_max; //For Richardsons Stationary
    double sigma = 0.75; //For SDslow
    double a_sds = 1; //For SDslow

    CJVector r_cg = qp.r(xk); //Initial Calculations for CG
    CJVector d_cg = r_cg;  //for CG
    double sig_cg = dot(r_cg,r_cg); //..for CG

    CJMatrix Hk(qp.n,qp.n); //For BFGS
    Hk.makeIdentity();

    start = chrono::system_clock::now();
    while((k<kmax) && (norm(qp.r(xk))>eps)){
        if (method == "Steepest Descent"){
            op_ctr += SD_Update(qp,xk);
        }
        else if(method == "Richardson's Stationary"){
            op_ctr += RS_Update(qp,xk,lambda_max);
        }
        else if(method == "SDslow"){
            op_ctr += SDS_Update(qp,xk,a_sds,sigma);
        }
        else if(method == "Conjugate Gradient"){
            op_ctr += CG_Update(qp,xk,r_cg,d_cg,sig_cg);
        }
        else if(method == "Gauss Southwell"){
            op_ctr += SD_Update(qp,xk);
        }
        else if(method == "Inexact Newton"){
            op_ctr += SD_Update(qp,xk);
        }
        else if(method == "BFGS"){
            op_ctr += BFGS_Update(qp,xk,Hk);
        }
        else{
            cout<<"Invalid Method"<<endl;
            return OptResult();
        }

        k++;
    }
    end = chrono::system_clock::now();
    ex_time = end-start;

    OptResult info;
    info.method = method;
    if(method == "SDslow")
        info.sigma = sigma;
    else if(method == "Richardson's Stationary")
        info.lambda_max = lambda_max;
    info.problem_type = "Quadratic s.p.d";
    info.dimension = qp.n;
    info.problem_ID = qp.name;
    info.num_iterations = k;
    info.num_ops = op_ctr;
    info.ex_time = ex_time;
    info.solution = xk;
    info.f_eval = qp.f(xk);
    return info;
}



void testResultStruct(){
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> ex_time;
    start = chrono::system_clock::now();
    end = chrono::system_clock::now();
    ex_time = end-start;
    OptResult test;
    test.method = "Steepest Descent";
    test.problem_type = "Quadratic s.p.d";
    test.dimension = 3;
    test.problem_ID = "#1, Convex";
    test.num_iterations = 100;
    test.num_ops = 1;
    test.ex_time = ex_time;
    cout<<test<<endl;
}



//Example #1, n=2, Quadratic
void ex_1_all(){
    size_t n = 2;
    CJMatrix G(n,n,SYMMETRIC);
    G(0,0)=16.0; G(0,1) = 0.0;
    G(1,0)=0.0;  G(1,1) = 4.0;
    CJVector c = {2.0, 2.0};
    CJVector x = {1.5, 1.5};
    double eps = 1e-9;

    QP_Problem q = {"Hello Optimization World",n,G,c};

    OptResult test1 = solve_QP(q, x,"Steepest Descent", eps);
    OptResult test2 = solve_QP(q, x, "Richardson's Stationary", eps);
    OptResult test3 = solve_QP(q, x, "SDslow", eps);
    OptResult test4 = solve_QP(q, x, "Conjugate Gradient", eps);
    OptResult test5 = solve_QP(q, x, "Gauss Southwell", eps);
    OptResult test6 = solve_QP(q, x, "BFGS", eps);
    cout<<test1<<endl;
    cout<<test2<<endl;
    cout<<test3<<endl;
    cout<<test4<<endl;
    cout<<test5<<endl;
    cout<<test6<<endl;
}

//Example #2, n=3, Quadratic

void ex_2_all(){
    size_t n = 3;
    CJMatrix G(n,n,SYMMETRIC);
    G(0,0)=16.0; G(0,1) = 0.0; G(0,2) = 0.0;
    G(1,0)=0.0;  G(1,1) = 4.0; G(1,2) = 0.0;
    G(2,0)=0.0;  G(2,1) = 0.0; G(2,2) = 6.0;
    CJVector c = {1, 2, 3};
    CJVector x = {1.5, 1.5, 1.5};
    double eps = 1e-9;

    QP_Problem q = {"3D Convex s.p.d Quadratic",n,G,c};

    OptResult test1 = solve_QP(q, x,"Steepest Descent", eps);
    OptResult test2 = solve_QP(q, x, "Richardson's Stationary", eps);
    OptResult test3 = solve_QP(q, x, "SDslow", eps);
    OptResult test4 = solve_QP(q, x, "Conjugate Gradient", eps);
    OptResult test5 = solve_QP(q, x, "Gauss Southwell", eps);
    cout<<test1<<endl;
    cout<<test2<<endl;
    cout<<test3<<endl;
    cout<<test4<<endl;
    cout<<test5<<endl;
}

void ex_3_all(){
    size_t n = 10;
    CJMatrix G(n,n,SYMMETRIC);
    CJVector c(n);
    for (int i =0; i<n;i++){
        G(i,i) = 5*sqrt(i);
        c[i] = 100/i;
    }
    CJVector x(n);



    double eps = 1e-9;

    QP_Problem q = {"n=10 Convex s.p.d Quadratic",n,G,c};

    OptResult test1 = solve_QP(q, x,"Steepest Descent", eps);
    OptResult test2 = solve_QP(q, x, "Richardson's Stationary", eps);
    OptResult test3 = solve_QP(q, x, "SDslow", eps);
    OptResult test4 = solve_QP(q, x, "Conjugate Gradient", eps);
    OptResult test5 = solve_QP(q, x, "Gauss Southwell", eps);
    cout<<test1<<endl;
    cout<<test2<<endl;
    cout<<test3<<endl;
    cout<<test4<<endl;
    cout<<test5<<endl;
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

    ex_1_all();
    ex_2_all();
    ex_3_all();
    //Experiments
        //Compare ALL Methods on a few 2-D, 3-D, & Large-Scale Problems

        //Examine effect of starting point X0 for ALL Methods

        //Examine effect of parameter SIGMA for SDSLOW

        //Examine effect of ALPHA for RICHARDSONS STATIONARY

    return 0;
}