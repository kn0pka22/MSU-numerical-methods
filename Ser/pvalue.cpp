#include "solution.hpp"
/*int p(int N){
    double *Nodes = new double[N+2];
    double *ValuesInNodes = new double[N+2];
    double *Phi = new double[N+1];
    double *AproximateValues = new double[N+2];
    double *Coef = new double[N+1];
    double *Nodes2 = new double[2*N+2];
    double *ValuesInNodes2 = new double[2*N+2];
    double *Phi2 = new double[2*N+1];
    double *Coef2 = new double[2*N+1];
    double *AproximateValues2 = new double[2*N+2];

    double h1 = 1/((double)N-1);
    double h2 = 1/(double)((2*N)-1);

    double h = h1;
    Nodes[0] = -h/2;
    Nodes[N+1] = 1 + h/2;
    for (int i = 1; i < N+1; i++){
        Nodes[i] = -h/2 + i*h;
    }
    for (int j = 0; j <= N+1; j++){
        ValuesInNodes[j] = f(Nodes[j]);
    }
    coefCalculating(Coef, Phi, ValuesInNodes, N);
    aproximateValuesCalculating(Coef, AproximateValues, Nodes, N);

    h = h2;
    Nodes2[0] = -h/2;
    Nodes2[2*N+1] = 1 + h/2;
    for (int i = 1; i < 2*N+1; i++){
        Nodes2[i] = -h/2 + i*h;
    }

    for (int j = 0; j <= 2*N+1; j++){
        ValuesInNodes2[j] = f(Nodes2[j]);
    }

    coefCalculating(Coef2, Phi2, ValuesInNodes2, 2*N);

    aproximateValuesCalculating(Coef2, AproximateValues2, Nodes2, 2*N);

    double err1 = normFunction(f, Coef, N);
    std::cout<<"err1 = "<<fabs(err1)<<std::endl;
    double err2 = normFunction(f, Coef2, 2*N);
    std::cout<<"err2 = "<<fabs(err1)<<std::endl;

    double a = log(err1/err2);
    double b = log(h1/h2);
    std::cout<<"p =  "<<fabs(a/b)<<std::endl;

    delete[] Nodes;
    delete[] ValuesInNodes;
    delete[] Phi;
    delete[] Coef;
    delete[] AproximateValues;
    delete[] Nodes2;
    delete[] ValuesInNodes2;
    delete[] Phi2;
    delete[] Coef2;
    delete[] AproximateValues2;
    return 0;
}
*/
    
    

