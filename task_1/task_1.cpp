#include "task_1.hpp"

double f(double x){
    return(2.*M_PI*x);
}


double scalar(double *x, double*y, unsigned int N){
    double res = 0.;
    for (unsigned int i=1;i<N;++i){
        res += x[i]*y[i];
    }
    return res;
}


int NodesGenerate(double *nodes, int N){
    double x_k = 0;
    double h = 0;

    h = 1 / ((double)N - 0.5);
    x_k = h/2;
    nodes[0] = - h / 2;
    nodes[N] = 1;

    for (int k = 1; k < N; k++){
        nodes[k] = x_k;
        x_k += h;
    }
    return 0;
}

double Phi(int m, int k, int h){
    return sin( M_PI * m * ( - h / 2 + k * h) );
}

void WritePhiTo(int m, int N, double* ph){
    double h = 1/(N - 0.5);
    for(int k = 1; k < N; ++k) {
        ph[k] = Phi(m, k, h);
    }
}


int CoeffCalculate(int N, double* c_m, double (*f)(double), double* x_k, double* u_k, double* phi){
    //double h = 1 / (N - 0.5);
    double c = 0;

    for(int i = 1; i < N; ++i){
        u_k[i] = f(x_k[i]);
    }

    for(int m = 1; m < N; ++m) {
        WritePhiTo(m, N, phi);
        c = 0;
        c = scalar(phi, u_k, N);// ans += phi[j] * u_k[j] * h;
        c_m[m] = 2 * c;
    }

    return 0;
}

double FourierCompute(double* coefs, int N, double x){
    double res = 0;
    for(int m = 1; m < N; ++m){
        res += coefs[m] * sin( M_PI * m * x );
    }
    return res;
}

double FullCompute(double x, int N, double* c_m, double (*f)(double), double* nodes, double* u_k, double* phi) {
    if(CoeffCalculate(N, c_m, f, nodes, u_k, phi) == 0){
        return FourierCompute(c_m, N, x);
    }
    std::cout<<"smth went wrong"<<std::endl;
    return 0;
}


