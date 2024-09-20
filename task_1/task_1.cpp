#include "task_1.hpp"


double f(double x){
    return(sin(M_PI * x));
}


void FillingNodes(double* xk, int N){

    double h = 1/(N-0.5);

    xk[0] = -h/2.;
    for (int i=1;i<N+1;++i){
        xk[i] = xk[i-1]+h;
    }

}

void FillingValues(double* xk, double* yk, double (*f)(double), int N){

    if (!xk) { 
        std::cerr<<"Please, fill nodes xk!"<<std::endl; 
        throw std::runtime_error("xk not found");
    }

    for (int i=0;i<N+1;++i){
        yk[i] = f(xk[i]);
    }

}

double PhiCalculate(int n, int k, double h){
    return sin(M_PI * n * (-h/2. + k*h));
}

void PhiVectorCalculate(int N, int n, double* phi){  //phi_k^{(n)} for fixed k
    double h = 1/(N-0.5);
    //phi[0] = 0;
    for (int k=1;k<N;++k){
        phi[k] = PhiCalculate(n, k, h);
    }
    //phi[N]=0;
}

double ScalarProduct(double* ar1,double* ar2, int N){
    double res = 0;
    for (int k=1;k<N; ++k){
        res+=ar1[k] * ar2[k];   //the original scalar product requires multiplication by h, but it cancels out
    }
    return res;
}

void CoeffCalculate(int N, int k, double* yk, double* phi, double* cn){
    double a = 0;
    double b = 0; 
    for (int n=1; n<N; ++n){
        PhiVectorCalculate(N, n, phi);
        a = ScalarProduct(yk, phi, N);
        b = ScalarProduct(phi,phi, N);
        cn[n] = a/b;
    }
}

double scalar(double *x, double*y, unsigned int N){
    double res = 0.;
    for (unsigned int i=1;i<N;++i){
        res += x[i]*y[i];
    }
    return res;
}

double FourierCompute(double* cn, int N, double x){
    double res = 0;
    for(int n = 1; n < N; ++n){
        res += cn[n] * sin( M_PI * n * x );
    }
    return res;
}
