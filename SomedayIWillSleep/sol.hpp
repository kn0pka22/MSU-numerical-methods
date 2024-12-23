#ifndef SOLUTION_H
#define SOLUTION_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>


# define M_PI 3.14159265358979323846
double f(double x);
double scalarProduct(double *Phi, double *ValuesInNodes, int N);
void phiCalculating(double *Phi, int m, int N);
double coefCalculating(double *Coef, double *Phi, double *ValuesInNodes, int N);
void aproximateValuesCalculating(double *Coef, double *AproximateValues, double *Nodes, int N);
void copyCoefToFile(int N, double *Coef);
void copyPointsToFile(int N, double *Nodes, double *ValuesInNodes, double *AproximateValues);
double backFourier(double *Coef, double x, int N);
double normFunction(double (*f)(double), double *Coef, int N);
int p(int N);
void printVector(double *Vector, int N);

double Lambda(int n, int N, double p);
void MatrixFill(double p, double* M, int N);
void printMatrix(double* matrix, int N);
void searchCoef(double *coef, double *b, double *phi, double p, int N);
void searchSol(double *coef, double *x, int N);
void MultiplyMatrixByVector(double* matrix, double* res, double* vec, int N);




#endif