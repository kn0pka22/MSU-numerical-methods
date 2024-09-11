#include "task_1.hpp"



double f(double x){
    return(2.*M_PI*x);
}

double scalar(double *x, double*y, unsigned int M){
    double res = 0.;
    double dx;
    for (unsigned int i=1;i<(M+1);++i){
        //dx = (x[i]+x[i-1])/2.;
        dx = x[i]-x[i-1];
        res += y[i]*dx;
    }

    return 0;
}