#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "linearSystemSolver.hpp"

#define MAX_ITERATIONS 15


double f(double x);
int GenerateEquidistantNodes(double a, double b, std::vector<double>& nodes);
int GenerateChebyshevNodes(double a, double b, std::vector<double>& nodes);
void FillingValues( std::vector<double>& nodes,  std::vector<double>& values, double (*f)(double), int N);
int MatrixFill(std::vector<double>& matrix, const std::vector<double>& nodes);
void printMatrix(const std::vector<double>& matrix);

void ExtendedNodes(std::vector<double>& ExNodes, std::vector<double>& Nodes);
void ExtendedValues(std::vector<double>& ExValues, std::vector<double>& Values);
void ExtendedF(std::vector<double>& ExF, const std::vector<double>& result, const std::vector<double>& ExNodes, int MM);
double delta(std::vector<double>& ExValues, std::vector<double>& ExF, int MM);
void CreateSigma(std::vector<double>& sigma, std::vector<double>& nodes, int MM, int N);

double MaxDeviation(const std::vector<double>& nodes, std::vector<double>& sigma, const std::vector<double>& coeffs, std::vector<double>& values, const std::vector<double>& valuesAll, int MM, int N);
double CalcPolynom(const std::vector<double>& result, double x, int N);
void  WriteToFile(double a, double b, const std::string& filename, std::vector<double>& coeffs);