#include "task_6.hpp"



int main(int argc, char* argv[]){

    int N;
    double xa, xb, ya, yb;

    if (argc<6 || argc>6){ std::cout<<"Please enter argc = 6!\n"; return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<1) || (sscanf(argv[2], "%lf", &xa)!=1) || (sscanf(argv[3], "%lf", &xb)!=1)
    || (sscanf(argv[2], "%lf", &ya)!=1) || (sscanf(argv[3], "%lf", &yb)!=1)){
        std::cout<<"Invalid input!\n \
        * N  – number of grid nodes, \n \
        * xa – left edge of the segment along the x-axis, \n \
        * xb – right edge of the segment along the x-axis, \n \
        * ya – left edge of the segment along the y-axis, \n \
        * yb – right edge of the segment along the y-axis, \n \
        Please enter:\n\
        N, xa, xb, ya, yb"<<std::endl;
        return -1;
    }
    //std::cout<<"N = "<<N<<" xa = "<<xa<<" xb = "<<xb<<" ya = "<<ya<<" yb = "<<yb<<std::endl;

    std::vector<FunctionWithName> functions = {
        {0, [](double x, double y) { return x + y; }, "f(x, y) = x + y"},
        {1, [](double x, double y) { return x * x + y * y; }, "f(x, y) = x^2 + y^2"},
        {2, [](double x, double y) { return x * x * x * x + x * x * y * y + y * y * y * y; }, "f(x, y) = x^4 + x^2 * y^2 + y^4"},
        {3, [](double x, double y) { return x - y; }, "f(x, y) = x - y"},
        {4, [](double x, double y) { return x * y; }, "f(x, y) = x * y"},
        {5, [](double x, double y) { return sin(x * y); }, "f(x, y) = sin(x * y)"}
    };
    int num=3;
    
    triangulation(N, xa, xb, ya, yb, "out.txt");
    std::cout<<"Quadrature =    "<<IntegrateQuadr1(N, xa, xb, ya, yb, functions[num].func, "out.txt")<<std::endl;
    std::cout<<"Analytical =    "<<functions[num].TrueRes(xa, xb, ya, yb)<<std::endl;
    
    int numTests = 6;
    GenereteFileForPCalculation(numTests, xa, xb, ya, yb, functions[num], N, "out.txt", "p.txt");
    
    return 0;
}

