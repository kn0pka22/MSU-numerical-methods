#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <functional>


double f(double x);
void printMatrix(const std::vector<double>& matrix);
void printVector(const std::vector<double>& vec);
std::vector<double> MultiplyMatrixByVector(const std::vector<double>& matrix, const std::vector<double>& vec);

void MatrixFill(double p, std::vector<double>& M);
void RightSideFill(const std::vector<double>& A, std::vector<double>& b); 

double psi(int k, int n, int N);
double Lambdan(int n, int N, double p);
double Dn(int n, std::vector<double>& f, double p, int N);
double FourierMethod(std::vector<double>& y, int N, double p, std::vector<double>& f);
double Richardson(std::vector<double>& x, const std::vector<double>& A, const std::vector<double>& b, double tau, int n, int mIter, std::vector<double>& mem);
double ErNorm(const std::vector<double>& A, const std::vector<double>& b, std::vector<double>& x, int N, std::vector<double>& mem);

