#include <iostream>
#include <fstream>
#include <cmath>
#include "solution.hpp"

void ClearFile(const std::string& filename) {
    std::ofstream ofs(filename, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
}

void WriteErrorToFile(int N, double error, const std::string& filename) {
    std::ofstream ofs(filename, std::ofstream::app);
    ofs << N << " " << error << std::endl;
}


int main() {
    int N = 10; 
    double p = 1.;

    const int numOfTests = 6;
    double errors[numOfTests]; 
    double h[numOfTests];

    ClearFile("error_data.txt");

    int NCopy = N;

    for (int test = 0; test < numOfTests; ++test) {
        double *pk = new double[N-1];
        for (int k = 0; k < N-1; k++) {
            pk[k] = 0.;
        }

        double *Nodes = new double[N+1];
        double *ValuesInNodes = new double[N+1];
        double *FullMatrix = new double[(N+1)*(N+1)];

        double *BasicNodes = new double[N-1];
        double *b = new double[N-1];
        double *B = new double[(N-1)*(N-1)];
        double *A = new double[(N-1)*(N-1)];

        double *Phi = new double[N-1];
        double *lambda = new double[N-1];
        double *Coef = new double[N-1];
        double *Ax = new double[N-1];
        double *memory1 = new double[N-1];
        double *x = new double[N-1];
        double *res = new double[N-1];
        double *res1 = new double[N-1];
        double *temp = new double[N-1];

        fillingNodes(Nodes, N);
        // std::cout << "---------------------------------------------------------------------------- Nodes ---------------------------------------------------------------------------" << std::endl;
        // printVector(Nodes, N+1);

        fillingValuesInNodes(Nodes, ValuesInNodes, N);
        // std::cout << "----------------------------------------------------------------------- Values in Nodes ----------------------------------------------------------------------" << std::endl;
        // printVector(ValuesInNodes, N+1);

        fillingBasicNodes(BasicNodes, N);
        // std::cout << "------------------------------------------------------------------------ Basic Nodes -------------------------------------------------------------------------" << std::endl;
        // printVector(BasicNodes, N-1);

        fillingValuesInBasicNodes(BasicNodes, b, N);
        // std::cout << "------------------------------------------------------------------------------ b -----------------------------------------------------------------------------" << std::endl;
        // printVector(b, N-1);

        fillingMatrixA(A, p, pk, N);
        // std::cout << "-------------------------------------------------------------------------- Matrix A --------------------------------------------------------------------------" << std::endl;
        // printMatrix(A, N-1);

        bool flag = searchCoef(Coef, b, Phi, p, N, lambda);
        if (!flag) {
            std::cout << "INVALID INPUT FOR SOLVING SYSTEM" << std::endl;
            return -1;
        }
        searchSol(Coef, x, N);

        errors[test] = ErNormInf(b, x, memory1, N, pk, p);
        h[test] = 1.0 / (N - 1.);
 
        if (test != 0) {
            if (std::fabs(h[test])>1e-10 && std::fabs(log(h[test - 1] / h[test]))>1e-10){
                double convergenceRate = log(errors[test - 1] / errors[test]) / log(h[test - 1] / h[test]);
                WriteErrorToFile(N, fabs(convergenceRate), "error_data.txt");
                // WriteErrorToFile(N,  errors[test], "error_data.txt");
                
            }
        }

        delete[] Nodes;
        delete[] ValuesInNodes;
        delete[] FullMatrix;
        delete[] BasicNodes;
        delete[] b;
        delete[] B;
        delete[] A;
        delete[] Phi;
        delete[] lambda;
        delete[] Coef;
        delete[] Ax;
        delete[] memory1;
        delete[] x;
        delete[] res;
        delete[] res1;
        delete[] temp;
        delete[] pk;

        N *= 2;
    }

    double h1 = 1.0 / (NCopy * pow(2, numOfTests - 2) - 1.);
    double h2 = 1.0 / (NCopy * pow(2, numOfTests - 1) - 1.);
    std::cout << "N1 = " << NCopy * pow(2, numOfTests - 2) << std::endl;
    std::cout << "N2 = " << NCopy * pow(2, numOfTests - 1) << std::endl;

    double a = log(errors[numOfTests - 2] / errors[numOfTests - 1]);
    double b = log(h1 / h2);
    double p_calculated = fabs(a / b);
    std::cout << " p =  " << p_calculated << std::endl;

    return 0;
}
