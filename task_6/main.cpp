#include "task_6.hpp"



int main(int argc, char* argv[]){
    int N;
    double Lx, Ly;
    if (argc<5 || argc>6){ std::cout<<"Please enter argc = 5 or 6!\n"; return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<2) || (sscanf(argv[2], "%lf", &Lx)!=1) || (sscanf(argv[3], "%lf", &Ly)!=1)){
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        Please enter:\n\
        N, Lx, Ly, k(>0) \n \
        or \n\
        N, Lx, Ly, 0 filename\n";
        return -1;
    }

    std::vector<FunctionWithName> functions = {
        {0, fun0, "f(x, y) = x+y"}
    };
    
    triangulation(N, Lx, Ly, "out.txt");
    std::cout<<"Quadrature =    "<<IntegrateQuadr1(N, fun1)<<std::endl;
    std::cout<<"Analytical =    "<<functions[0].TrueRes(0, Lx, 0, Ly)<<std::endl;

    
    return 0;
}

