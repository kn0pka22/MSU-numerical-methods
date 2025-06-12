#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
// # define M_PI          3.14159265358979323846


double f(double t, double x, double y);
// double f(double t, double x, double y, double p);

double u0(double x, double y);
double uFunc(double t, double x, double y);

double u0WithT(double t, double x, double y);
void CalculateFourierCoefficients(int Nx, int Ny, double t, double* xk, double* yk, 
                                double* level, double* D, double* C, 
                                double* phi, double* fMem, 
                                double (*f)(double, double, double));

void FillingNodes(double* xk, int N);
// void PrintMatrix(const double* matrix, int N, const std::string& name);
void PrintMatrix(const double* matrix, int rows, int columns, const std::string& name);
//void FillingValues(double* xk, double* yk, double (*f)(double), int N);
void FillingUMatrix(int N, double* U, double* xk, double (*f)(double, double));
double PhiCalculate(int n, int k, double h);
void   PhiVectorCalculate(int N, int k, double* phi);
double ScalarProduct(double* ar1,double* ar2, int N);
void   CoeffCalculate(int N, double* yk, double* phi, double* cn);
double FourierCompute(double* cn, int N, double x);

// void FillingDMatrix(int N, double* D, double* U, double* phi);  
void FillingDMatrix(int Nx, int Ny, double* D, double* U, double* phi);
void FillingCMatrix(int Nx, int Ny, double* D, double* C, double* fMem, double* phi);
// void FillingCMatrix(int N, double* D, double* C, double *fmemory, double* phi);

// double Calc2DFourier(double* C, int N, double x, double y);
double Calc2DFourier(double* C, int Nx, int Ny, double x, double y);

double WriteToConsole(int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi);

// void  WriteToFile(const std::string& filename, int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi);
// void  WriteToFile(std::string flag, const std::string& filename, int Nx, int Ny, int Nt, double* tk, double* xk, double* yk, double* U, double* C, double* D, double* fmemory, double* phi);

// void WriteToFile(std::string flag, const std::string& filename, int Nx, int Ny, int Nt, 
//                 double* tk, double* xk, double* yk, double* U, double* C, double* D, 
//                 double* fmemory, double* phi, bool onlyNodes = 1);

//  void WriteToFile(std::string flag, const std::string& filename, int Nx, int Ny, int Nt, 
//                 double* tk, double* xk, double* yk, double* U, double* C, double* D, 
//                 double* fmemory, double* phi, 
//                 double (*f_func)(double, double, double), bool onlyNodes = 1);

void WriteToFileSimple(const std::string& filename, int Nx, int Ny, 
                double tFinal, double* xk, double* yk, double* C, double* D, 
                double* fmemory, double* phi, 
                double (*f_func)(double, double, double)) ;


int pcalculate(int N);
// double normFunction(double (*f)(double, double), double *C, int N);


// void coeffsForf(double (*f)(double, double), double* U, double* xk, int N, double *C, double tau, int layer, double* phi, double* D, double* fmemory, double* d_mn);
// void FillingU0Matrix(int N, double* U, double* xk, double t, double (*f)(double, double));
// void FillingU0Matrix(int N, double* U, double* xk, double (*u0)(double, double));
void FillingU0Matrix(int Nx, int Ny, double* U, double* xk, double* yk, double (*u0)(double, double));
// void FillingUMatrix(int N, double* U, double* xk, double (*f)(double, double), double* U_layer);

void FillingFMatrix(int Nx, int Ny, double* F, double t, double* xk, double* yk, double (*f)(double, double, double));
// void FillingFMatrix(int N, double* F, double t, double* xk, double* yk, double (*f)(double, double, double));
void FillingNodesForTime(double* tk, int Nt);
// void FillingUMatrix(int N, double* U, double* xk, double (*u)(double, double));

void FillingUMatrix(int N, double* U, double* xk, double t, double (*u)(double, double, double));
double Calc1DFourier(double* C, int N, double x);


// void FindFourierCoefs(double* U, double* D, double* C, int N,
//                       double* fMem, double* xk, double (*u)(double, double),
//                       double* netmemory, double* uMem, double* phi);

// void solvePDE2d(double (*u0)(double, double),
//                 double (*f)(double, double, double), 
//                 int N, int Nt,
//                 double* uijn,
//                 double* U, double* D, double* Cmatrix,
//                 double* xk, double* fMem, double* phi);                      

double lambda(int N, double n);
// void CalculateCoeffsForNextLevel(int Nx, int Ny, double tau, double* fCoeffs, double* resCoeffs);
void CalculateCoeffsForNextLevel(int Nx, int Ny, double p, double tau, double* fCoeffs, double* resCoeffs);
void AMultiplyByX(int Nx, int Ny, double tau, double* xk, double* yk, double* u, double* res);

void CalculateFourierCoefficientsWithReadyMatrix(int Nx, double Ny, double* matr, double* D, double* C, 
                                            double* phi, double* fMem);
double error(int Nx, int Ny, int Nt, double* xk, double* yk, double* u, double* b, double* prod);

void BSolver(int Nx, int Ny, int Nt, double p, double theta, 
            double* u, double* xk, double* yk, double* b, double* coeffsForU, double* coeffsForF,
            double* prod, double* y, double* phi, double* fMem, double* D);

void WriteToFileSimple(const std::string& filename, int Nx, int Ny, 
                double tFinal, double* xk, double* yk, double* U, double* C, double* D, 
                double* fmemory, double* phi, 
                double (*f_func)(double, double, double));

void WriteSolutionWithError(int Nx, int Ny, double t, double* xk, double* yk, 
                          double* level, const std::string& filename);

double uuuFunc(double* b, int i, int j, int Nx);