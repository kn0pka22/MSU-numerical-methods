#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
//#include <regex>


struct FunctionWithName{
    int id;
    std::function<double(double)> func;
    std::string name;
    double a;
    double b;
    double TrueRes();    
};

//double f(double x);
double f0(double x);
double f1(double x);
double f2(double x);
double f3(double x);
double f4(double x);
double f5(double x);
double f6(double x);
double f7(double x);
double f8(double x);
double f9(double x);

double g0(double x);
double g1(double x);
double g2(double x);
double g3(double x);

double h0(double x);
double h1(double x);
double h2(double x);
double h3(double x);
double h4(double x);
double h5(double x);
double h6(double x);
double Runge(double x);

double f14(double x);
double f15(double x);
double f16(double x);
double f17(double x);
double f18(double x);
double f19(double x);
double f20(double x);
double f21(double x);
double f22(double x);
double f23(double x);
double f24(double x);
double f25(double x);
double f26(double x);
double f27(double x);
double f28(double x);
double f29(double x);
double f30(double x);


double RectangleMethod(double a, double b, const std::function<double(double)>& f);

double SimpsonMethod(double a, double b, const std::function<double(double)>& f);

double GaussMethod(double a, double b, const std::function<double(double)>& f);

double CompositeRectangleMethod(double a, double b,const std::function<double(double)>& f, int N);

double CompositeSimpsonQuadrature (double a, double b, const std::function<double(double)>& f, int N);

double CompositeGaussianQuadrature(double a, double b, const std::function<double(double)>& f, int N);

void writeResultsToFile(const std::string& filename, const std::vector<FunctionWithName>& functions, std::vector<double>& WolframResults);
void writeResultsToFile(const std::string& filename, std::vector<FunctionWithName>& functions, int N);
void GenereteFileForPCalculation(int numTests, int numFunc, int N, const std::string& filename, std::vector<FunctionWithName>& functions);