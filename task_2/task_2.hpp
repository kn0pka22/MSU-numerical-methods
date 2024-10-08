#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
# define M_PI          3.14159265358979323846


double f(double x, double y);
void FillingNodes(double* xk, int N);
void PrintMatrix(const double* matrix, int N, const std::string& name);
//void FillingValues(double* xk, double* yk, double (*f)(double), int N);
void FillingUMatrix(int N, double* U, double* xk, double (*f)(double, double));
double PhiCalculate(int n, int k, double h);
void   PhiVectorCalculate(int N, int k, double* phi);
double ScalarProduct(double* ar1,double* ar2, int N);
void   CoeffCalculate(int N, double* yk, double* phi, double* cn);
double FourierCompute(double* cn, int N, double x);

void FillingDMatrix(int N, double* D, double* U, double* phi);  
void FillingCMatrix(int N, double* D, double* C, double *fmemory, double* phi);

double Calc2DFourier(double* C, int N, double x, double y);

void WriteToConsole(int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi);
void  WriteToFile(const std::string& filename, int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi);
int pcalculate(int N);
