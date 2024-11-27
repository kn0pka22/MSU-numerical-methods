#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <string>

struct FunctionWithName{
    int id;
    std::function<double(double, double)> func;
    std::string name;
    double TrueRes(double xa, double xb, double ya, double yb);    
};


double fun0(double x, double y);
double fun1(double x, double y);
double fun2(double x, double y);
double fun3(double x, double y);
double fun4(double x, double y);

void triangulation(int N, double Lx, double Ly, const std::string &filename);
double IntegrateQuadr1(int N, std::function<double(double, double)> f);
double IntegrateQuadr2(int N, std::function<double(double, double)> f);