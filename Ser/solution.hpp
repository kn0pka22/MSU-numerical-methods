#ifndef SOLUTION_H
#define SOLUTION_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>


# define M_PI 3.14159265358979323846


double f(double x, double y);
void fillingNodes(double *Nodes, int N);
void PrintMatrix(const double* matrix, int N, const std::string& name);
void fillingValuesInNodes(double *Nodes, double *ValuesInNodes, int N);
void printNodes(double *Nodes, int N);
void printMatrix(double *matrix, int N);
void printValueInNodes(double *ValuesInNodes, int N);
double scalarProduct(double *Phi, double *ValuesInNodes, int N);
void phiCalculating(double *Phi, int m, int N);
void coefCalculating(double *Coef, double *Phi, double *ValuesInNodes, int N);
//void fillingCMatrix();
void fillingCMatrix(double *C, double *Phi, double *memory, double *D, int N);
void fillingDMatrix(double* D, double* Phi, double* ValuesInNodes, int N);
double WriteToConsole(int N, double* xk,  double* ValuesInNodes, double* C, double* D, double* Memory, double* Phi);
void aproximateValuesCalculating(double *Coef, double *AproximateValues, double *Nodes, int N);
void copyCoefToFile(int N, double *Coef);
void copyPointsToFile(int N, double *Nodes, double *ValuesInNodes, double *AproximateValues);
double backFourier(double *C, double x, double y, int N);
double normFunction(double (*f)(double), double *Coef, int N);
//int p(int N);

#endif