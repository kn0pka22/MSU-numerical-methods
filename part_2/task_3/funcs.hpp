#ifndef FUNCS_HPP
#define FUNCS_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>


double func(double t, double x);
double u0(double x);
 
void solveImplicit(int M, double levels, double tau, double** u);

int findSteps(double T,double tau);
void createThreeDiag(double* upDiag, double* midDiag, double* downDiag, int M, double h, double tau);
void SolveTridiagonal(double* upDiag, double* midDiag, double* downDiag, double* rhs, int M);   // метод прогонки                 
void freeSolution(double** u, int M);
void printMatr(double** matr, int N, int M);
void printVec(double* vec, int N);
void solveExplicit(int N, int M, double tau, double h, double** u);

double error(double** u, int M, int N, double h, double T);
void runConvergenceTest(const std::string& filename, double T);

void saveMatrixToFile(double** matr, int rows, int cols, const std::string& filename);
void saveVectorToFile(double* vec, int size, const std::string& filename);


#endif