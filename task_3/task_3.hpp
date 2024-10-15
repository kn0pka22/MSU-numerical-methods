#include "linearSystemSolver.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

# define M_PI           3.14159265358979323846


double f(double x);
int GenerateEquidistantNodes(double a, double b, double (*f)(double), std::vector<double>& knots, std::vector<double>& values);
int GenerateChebyshevNodes(double a, double b, double (*f)(double), std::vector<double>& knots, std::vector<double>& values);
int MatrixFill(std::vector<double>& matrix, const std::vector<double>& nodes);
double PnCalculation(std::vector<double>& coeffs, double x);
double LnCalculation(const std::vector<double>& x, const std::vector<double>& y, double x_value);
int VecFill(std::vector<double>& vec, const std::vector<double>& values);
void  WriteToFilePn(double a, double b, const std::string& filename, std::vector<double>& coeffs);
void  WriteToFileLn(double a, double b, const std::string& filename, std::vector<double>& nodes, std::vector<double>& values);




