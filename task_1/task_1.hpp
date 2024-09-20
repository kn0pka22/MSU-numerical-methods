#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>




double f(double x);
double scalar(double *x, double*y, unsigned int N);
void   FillingNodes(double* xk, int N); //генерация узлов
void   FillingValues(double* xk, double* yk, double (*f)(double), int N);
double PhiCalculate(int n, int k, double h);
void   PhiVectorCalculate(int N, int k, double* phi);
double ScalarProduct(double* ar1,double* ar2, int N);
void   CoeffCalculate(int N, int k, double* yk, double* phi, double* cn);
double FourierCompute(double* cn, int N, double x);


