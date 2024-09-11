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


int GrigGenerate(double *grid, int N) 
{
    double x_k = 0;
    double h = 0;

    if (N < 3)
        return -1;
    h = 1 / ((double)N - 0.5);
    x_k = h / 2;
    grid[0] = - h / 2;
    grid[N] = 1;

    for (int k = 1; k < N; k++)
    {
        grid[k] = x_k;
        x_k += h;
    }
    return 0;
}