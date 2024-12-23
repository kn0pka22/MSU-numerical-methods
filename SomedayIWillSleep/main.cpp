#include<iostream>
#include <fstream>
#include "sol.hpp"

double f(double x){
    return cos(M_PI * 7. * x);
    //return 1 + 0. * x;
    //return 2./3. * x*x*x - x*x;
    //return 2./3. * pow(x, 12) - pow(x, 8);
    //return x*x * (1-x)*(1-x) * exp(-x);
}

int main(){

    int N = 5;
    double p=0.;
    double h = 1. / (double)(N-1);

    double *Nodes = new double[N+1];
    double *ValuesInNodes = new double[N+1];
    double *Phi = new double[N+1];
    double *Coef = new double[N+1];
    double *AproximateValues = new double[N+1];
    double *Matrix = new double[(N+1)*(N+1)];
    double *res = new double[N+1];



    Nodes[0] = -h/2.;
    Nodes[N] = 1 + h/2.;
    for (int i = 1; i < N+1; i++){
        Nodes[i] = -h/2. + i*h;
    }

    for (int j = 0; j < N+1; j++){
        ValuesInNodes[j] = f(Nodes[j]);
    }



    //coefCalculating(Coef, Phi, ValuesInNodes, N);

    //aproximateValuesCalculating(Coef, AproximateValues, Nodes, N);

    searchCoef(Coef, ValuesInNodes, Phi, p, N);
    searchSol(Coef, AproximateValues, N);

    MatrixFill(p, Matrix, N);
    MultiplyMatrixByVector(Matrix, res, AproximateValues, N+1);

    //====================CHECK======================================
    std::cout<<"CHECK: "<<std::endl;  // <- correct
    printVector( ValuesInNodes, N+1);
    printVector( res, N+1);
    //===============================================================


    //printMatrix(Matrix, N);

    


    //====================CHECK for 1d======================================
    // std::cout<<"CHECK: "<<std::endl;  // <- correct
    // std::cout<<"backFourier(Coef, 0.1, N) = " \
    //          <<backFourier(Coef, 0.1, N)<<std::endl;

    // std::cout<<"        f(0.1)            = "  \
    //          <<f(0.1)<<std::endl;
    //===============================================================




    delete[] Nodes;
    delete[] ValuesInNodes;
    delete[] Phi;
    delete[] Coef;
    delete[] AproximateValues;


    delete[] Matrix;
    delete[] res;
    return 1;
}
