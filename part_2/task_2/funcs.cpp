#include "funcs.hpp"

double ans(double A, double y0, double x){
    double trueSol = y0* exp(- A * y0 * x);
    return trueSol;
}

double scheme_1(double y0, double A, int N){ 
    double h = 1./N;
    double y, yNext;
    y = y0;
    double f, sum = 0.;

    for(int k = 1; k < N; k++){
        yNext = y * (1 - A*h);
        y = yNext;

        f = ans(A, y0, k * h);
        sum = sum + (f - y) * (f - y);
    }

    sum = sum*h;
    sum = sqrt(sum);

    return sum; 
}

double scheme_2(double y0, double A, int N){
    double h = 1./N;
    double y, yNext;
    y = y0;
    double f, sum = 0.;

    for(int k = 1; k < N; k++){
        yNext = y / (1. + A * h);
        y = yNext;

        f = ans(A, y0, k* h);
        sum = sum + (f - y) * (f - y);
    }
    sum = sum*h;
    sum = sqrt(sum);

    return sum; 
}

double scheme_3(double y0, double A, int N){
    double h = 1./N;
    double y, yNext;
    y = y0;
    double f, sum = 0.;

    for(int k = 1; k < N; k++){
        yNext = y * (2. - A * h) / (2 + A * h);
        y = yNext;

        f = ans(A, y0, k* h);
        sum = sum + (f - y) * (f - y);
    }
    sum = sum*h;
    sum = sqrt(sum);

    return sum; 
}

double scheme_4(double y0, double A, int N){
    double h = 0.01 / N;
    double y_old, y, yNext;
    y_old = y0;
    y = y0 * (1 - A * h);
    double f, sum = 0.;

    for(int k = 2; k < N; k++){
        yNext = y_old - 2. * A * h * y;
        y_old = y;
        y = yNext;

        f = ans(A, y0, k* h);
        sum = sum + (f - y) * (f - y);
    }
    sum = sum*h;
    sum = sqrt(sum);

    return sum; 
}

