#include<iostream>
#include <fstream>
#include "solution.hpp"

int main(){

    int N = 5;
    // double h = 1. / (double)(N-1);
    double p = 1.;
    // int mIter = 300; 

    double *Nodes = new double[N+1];
    double *ValuesInNodes = new double[N+1];
    double *FullMatrix = new double[(N+1)*(N+1)];

    double *BasicNodes = new double[N-1];
    double *ValuesInBasicNodes = new double[N-1];
    double *BasicMatrix = new double[(N-1)*(N-1)];

    double *Phi = new double[N-1];

    double *lambda = new double[N-1];
    
    double *Coef = new double[N-1];

    double *pk = new double[N-1];

    double *Ax = new double[N-1];
    double *memory1 = new double[N-1];

    double *x = new double[N-1];
    double *res = new double[N-1];
    double *res1 = new double[N-1];

    double *temp = new double[N-1];

    fillingNodes(Nodes, N);
    // std::cout << "Nodes:" << std::endl;
    // printVector(Nodes, N+1);
    // std::cout<<std::endl;

    fillingValuesInNodes(Nodes, ValuesInNodes, N);
    // std::cout << "Values in Nodes:" << std::endl;
    // printVector(ValuesInNodes, N+1);
    // std::cout<<std::endl;

    fillingBasicNodes(BasicNodes, N);
    // std::cout << "Basic Nodes:" << std::endl;
    // printVector(BasicNodes, N-1);
    // std::cout<<std::endl;

    fillingValuesInBasicNodes(BasicNodes, ValuesInBasicNodes, N);
    // std::cout << "Values in Basic Nodes:" << std::endl;
    // printVector(ValuesInBasicNodes, N-1);
    // std::cout<<std::endl;

    std::cout << "----------------------Basic Matrix----------------------------------------------" << std::endl;
    fillingBasicMatrix(p, BasicMatrix, N);
    printMatrix(BasicMatrix, N-1);

    std::cout << "-----------------------------Full Matrix-----------------------------------------" << std::endl;
    fillingFullMatrix(p, FullMatrix, N);
    printMatrix(FullMatrix, N+1);
    std::cout << "---------------------------------------------------------------------" << std::endl;

    Lambda(lambda, p, N);
 
    bool flag = searchCoef(Coef, ValuesInBasicNodes, Phi, p, N, lambda);
    if (!flag){
        std::cout<<"INVALID INPUT FOR SOLVING SYSTEM"<<std::endl;
        return -1;
    }
    searchSol(Coef, x, N); 

    for (int k = 0; k < N-1; k++){
        pk[k]=0.;
    }
    MultiplicationByA(x, temp, N, pk, p);


    //====================CHECK======================================
    std::cout<<"CHECK: "<<std::endl;  
    printVector(ValuesInBasicNodes, N-1);
    printVector(temp, N-1);
    //===============================================================

    std::cout<<"\n !!!!GENIOUS!!!! \n"<<std::endl;

    //====================================================================================
    int numberTest = 200;

    double mean = 0;
    for (int k = 0; k < N-1; k++){
        pk[k] = 0.;
        mean += pk[k] + p;
    }
    mean /= (double)(N-1.);

    Lambda(lambda, p, N);


    double EigenValueMin = lambda[0]; // Наименьшее с.з.
    double EigenValueMax = lambda[N-2]; // Наибольшее с.з.
    
    //double tau = 2. / (EigenValueMin + EigenValueMax);

    double tau = 1;

    double q = (EigenValueMax - EigenValueMin) / (EigenValueMax + EigenValueMin);
    //double q = SearchQ(BasicMatrix, N);
    double dq = 1;

    double q0 = q;
    double resid0 = BSolver(x, ValuesInBasicNodes, tau, 
                            1, res, res1, p, N, pk, Coef, Phi, 
                            lambda); //for x^0
    
    
    //std::cout << "resid0 = " << resid0 << std::endl;
    std::ofstream fout("out.txt");

    fout << "1  " << resid0 << " " << q * resid0 << "\n";
    double resid;

    q = q0;    
    std::cout << "q0 = " << q0 << std::endl;
    
    fout << std::fixed << std::setprecision(15);
    for (int iteration = 2; iteration < numberTest+2; iteration++){
        resid = BSolver(x, ValuesInBasicNodes, tau, iteration,
        res, res1, p, N, pk, Coef, Phi, lambda); 
        //std::cout << "for iter = " << iter << " started" << std::endl;
        fout << iteration+1 << " " << resid << " " << q0 * resid0 << "\n";
        q0 *= q0;
    }
    fout.close();



    
    delete[] Nodes;
    delete[] ValuesInNodes;
    delete[] FullMatrix;

    delete[] BasicNodes;
    delete[] ValuesInBasicNodes;
    delete[] BasicMatrix;

    delete[] Coef;
    delete[] Phi;
    delete[] lambda;

    delete[] res;
    delete[] res1;

    delete[] pk;

    delete[] temp;

    delete[] Ax;
    delete[] memory1;

    delete[] x;

    return 1;
}
