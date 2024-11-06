#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "linearSystemSolver.hpp"



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


double MaxDeviation(const std::vector<double>& nodes, std::vector<double>& sigma, const std::vector<double>& res, const std::vector<double>& valuesAll, int MM, int N);
//double CalcP(std::vector<double>& coeffs,  double x);