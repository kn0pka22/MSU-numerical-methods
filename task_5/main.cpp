#include "task_5.hpp"



int main(int argc, char* argv[]){
    int N;
    if (argc<2 || argc>2){ std::cout<<"Please enter argc = 2!\n"; return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<2)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        Please enter:\n\
        N, a, b, k\n";
        return -1;
    }

    double a = -1.;
    double b = 1.;

    std::vector<FunctionWithName> functions = {
        {0, f0, "f(x) = 1", a, b},
        {1, f1, "f(x) = 2x", a, b},
        {2, f2, "f(x) = 3x^2", a, b},
        {3, f3, "f(x) = 4x^3", a, b},
        {4, f4, "f(x) = e^x", 0., 1.},
        {5, f5, "f(x) = x * e^x", 0., 2.},
        {6, f6, "f(x) = cos(100x)", -M_PI, M_PI}, 
        {7, f7, "f(x) = e^(-10x)", 0., 1.},
        {8, f8, "f(x) = 1/sqrt(x)", 0.+1e-2, 1.},
        {9, f9, "f(x) = 1/sqrt(1 - x^2)", -1.+1e-4, 1.-1e-4},
        {10, g0, "g(x) = -10", a, b},
        {11, g0, "g(x) = -22x", a, b},
        {12, g0, "g(x) = -36x^2", a, b},
        {13, g0, "g(x) = -52x^3", a, b},
        {14, f14, "f(x) = sin(x)", -M_PI, 2*M_PI},
        {15, f15, "f(x) = cos(x)", -M_PI, M_PI},
        {16, f16, "f(x) = tan(x)", -M_PI/2 + 1e-4, M_PI/2 - 1e-4}, 
        {17, f17, "f(x) = log(x)", 0.1, 10.},  
        {18, f18, "f(x) = x^2 - 4x + 3", -5., 5.},  
        {19, f19, "f(x) = x^3 - 3x^2 + 2x", -2., 2.},  
        {20, f20, "f(x) = x * sin(x)", -M_PI, M_PI}, 
        {21, f21, "f(x) = e^(-x^2)", -5., 5.},  
        {22, f22, "f(x) = sqrt(x)", 0., 10.},  
        {23, f23, "f(x) = 1/x", 0.1, 10.},  
        {24, f24, "f(x) = exp(x) - 1", -1., 1.}, 
        {25, f25, "f(x) = sin(x)/x", 1e-4, 10.}, 
        {26, f26, "f(x) = x * cos(x)", -M_PI, M_PI},  
        {27, f27, "f(x) = x^4", -2., 2.},  
        {28, f28, "f(x) = sqrt(1 - x^2)", -1., 1.},  
        {29, f29, "f(x) = log(1 + x)", 0., 5.},  
        {30, f30, "f(x) = x^2 * e^x", 0., 2.},
        {31, Runge, "f(x) = 1/(1 + 25x^2)", -1.,1.}
};


    // std::ifstream inputFile("integralsFromWolfram.txt"); 
    // if (!inputFile) {
    //     std::cerr << "Error opening file!" << std::endl;
    //     return 1;
    // }

    // std::vector<double> WolframResults;
    // double WolframValue;
    // while (inputFile >> WolframValue) {
    //     WolframResults.push_back(WolframValue);
    // }
    // inputFile.close();
    //writeResultsToFile("results.txt", functions, WolframResults);


    writeResultsToFile("results.txt", functions, N);
    GenereteFileForPCalculation(10, 31, N, "output.txt", functions);

    
    return 0;
}

