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
        {0, f0, "f(x) = 1", a,b},
        {1, f1, "f(x) = 2x", a,b},
        {2, f2, "f(x) = 3x^2", a,b},
        {3, f3, "f(x) = 4x^3", a,b},
        {4, f4, "f(x) = e^x", 0.,1.},
        {5, f5, "f(x) = x * e^x", 0.,2.},
        {6, f6, "f(x) = cos(100x)", -M_PI,M_PI}, 
        {7, f7, "f(x) = e^(-10x)", 0., 1,},
        {8, f8, "f(x) = 1/sqrt(x)", 0.+1e-4, 1.},
        {9, f9, "f(x) = 1/sqrt(1 - x^2)", -1.+1e-4, 1.-1e-4},
        {10,g0, "g(x) = -10", a,b},
        {11,g0, "g(x) = -22x", a,b},
        {12,g0, "g(x) = -36x^2", a,b},
        {13,g0, "g(x) = -52x^3", a,b},
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


    writeResultsToFile("results.txt", functions);
    GenereteFileForPCalculation(10, 5, N, "output.txt", functions);

    
    return 0;
}

