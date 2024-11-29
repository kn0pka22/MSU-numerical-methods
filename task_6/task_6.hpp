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

double PointFromNumVert(int numVert, int N, double h, char c);
struct Triangle{
    //we number the vertices of the triangle counterclockwise, starting from the lower left vertex (if there is none, then from the lower right vertex)
    // . _ . \ . |
    // . | . - . "\"
    int v1, v2, v3;  
    Triangle(int v1, int v2, int v3);  
};

 
void triangulation(int N, double xa, double xb, double ya, double yb, const std::string &filename);
double IntegrateQuadr1(int N, double xa, double xb, double ya, double yb, std::function<double(double, double)>& f, const std::string& filename);
double IntegrateQuadr2(int N, double xa, double xb, double ya, double yb, std::function<double(double, double)>& f, const std::string& filename);

double GenereteFileForPCalculation(int numTests, double xa, double xb, double ya, double yb, 
                                FunctionWithName f, int N,
                                const std::string& fileForTriangulation,
                                const std::string& fileForP);