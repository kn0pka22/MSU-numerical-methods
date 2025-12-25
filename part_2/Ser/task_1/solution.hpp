#ifndef SOLUTION_H
#define SOLUTION_H
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>

// # define M_PI 3.14159265358979323846

double f(double x);

void fillingNodes(double *Nodes, int N);
void fillingValuesInNodes(double *Nodes, double *ValuesInNodes, int N);
void fillingFullMatrix(double p, double* M, int N);

void fillingBasicNodes(double *BasicNodes, int N);
void fillingValuesInBasicNodes(double *BasicNodes, double *ValuesInBasicNodes, int N);
void fillingMatrixB(double* B, double p, int N);
void fillingMatrixA(double *A, double p, double *pk, int N);

void MultiplicationByB(double *vector, double *res, int N, double p);
void MultiplicationByA(double *vector, double *res, int N, double p, double *pk);
// void MultbyNum(double *vector, double *res, double *lambda, int N, int m);

void Lambda(double *lambda, double p, int N);
void phiCalculating(double *Phi, int m, int N);
double scalarProduct(double *Phi, double *Vector, int N);
double coefCalculating(double *Coef, double *Phi, double *Vector, int N);


void printMatrix(double *Matrix, int N);
void printVector(double *Vector, int N);

bool searchCoef(double *coef, double *b, double *phi, double p, int N, double *lambda);
void searchSol(double *coef, double *x, int N);
double BSolver( double* x, double* b, 
                double tau, int mIter, double* mem, 
                double* mem1, double p, int N, double *pk,
                double *Coef, double *Phi, double *lambda);
double ErNormInf(double* b, double* x, double* mem, int N, double *pk, double p);
#endif