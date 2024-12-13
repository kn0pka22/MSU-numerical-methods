#include "task_1.hpp"

int pcalculate(int N){

double h;

double* xk  = new double[N+1];
double* xk2 = new double[2*N+1];

double* yk  = new double[N+1];
double* yk2   = new double[2*N+1];


double* phi  = new double[N+1];
double* phi2  = new double[2*N+1];

double* cn  = new double[N+1];
double* cn2 = new double[2*N+1];

FillingNodes(xk, N);
FillingValues(xk, yk, f, N);

FillingNodes(xk2, 2*N);
FillingValues(xk2, yk2, f, 2*N);

h = 1/(N-0.5);
for (int i = 1; i < N+1; ++i){
    CoeffCalculate(N, i, yk, phi, cn);
    FourierCompute(cn, N, (- h/2. + i* h));
}

h = 1/(2*N-0.5);
for (int i = 1; i < 2*N+1; ++i){
    CoeffCalculate(2*N, i, yk2, phi2, cn2);
    FourierCompute(cn2, 2*N, (- h/2. + i* h));
}

double err1 = NormFunction(f, cn,   N);
double err2 = NormFunction(f, cn2, 2*N);

double h1 = 1/((double)N-0.5);
double h2 = 1/((double)(2*N)-0.5);

double a = log(err1/err2);
double b = log(h1/h2);



//std::cout<<"(h1/h2)^p = err1/err2:"<<"  ";
std::cout<<"p =  "<<fabs(a/b)<<std::endl;


delete[] cn; 
delete[] cn2;
delete[] xk; 
delete[] xk2; 
delete[] yk;
delete[] yk2; 
delete[] phi; 
delete[] phi2; 

return 0;

}