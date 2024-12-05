#include "task_7.hpp"


int main(int argc, char *argv[]){
    int N;
    double p;

    if (argc<3 || argc>3){ std::cout<<"Please enter argc = 3!\n"; return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<1) || (sscanf(argv[2], "%lf", &p)!=1) || (p<0)){
        std::cout<<"Invalid input!\n \
        * N   – number of grid nodes, \n \
        * p>0 – parameter\n \
        Please enter:\n\
        N, p"<<std::endl;
        return -1;
    }
    

    std::vector<double> Nodes(N + 1);
    std::vector<double> BasicNodes(N - 1);
    
    std::vector<double> ValuesInNodes(N + 1);
    std::vector<double> ValuesInBasicNodes(N - 1);

    std::vector<double> BasicMatrix((N - 1) * (N - 1));
    std::vector<double> FullMatrix((N + 1) * (N + 1));

    std::vector<double> x(N - 1);
    std::vector<double> Fullx(N + 1);

    std::vector<double> mem(N - 1 );





    FillingNodes(Nodes);
    FillingBasicNodes(BasicNodes);
    FillingValues(Nodes, ValuesInNodes, f);
    FillingValues(BasicNodes, ValuesInBasicNodes, f);
    //printVector(ValuesInBasicNodes);
 
    BasisMatrixFill(p, BasicMatrix);
    FullMatrixFill(p, FullMatrix);
    //printMatrix(FullMatrix);
    //printMatrix(Matrix);


    Fourier(x, p, ValuesInBasicNodes);
    Fullx[0] = -x[0];
    for (int i=0;i<N-1;++i){
        Fullx[i+1] = x[i];
    }
    Fullx[N] = 0.;


 
    //std::cout<<"check: "<<std::endl; <- correct
    //printVector(MultiplyMatrixByVector(FullMatrix, Fullx));
    //printVector(ValuesInNodes);
 
//============================================================
    
    double EigenValueMin = Lambda(1, N, p);
    double EigenValueMax = Lambda(N-1, N, p);
    
    double tau = 2. / (EigenValueMin + EigenValueMax);
    double q = (EigenValueMax - EigenValueMin) / (EigenValueMax + EigenValueMin);
    double dq = 1;


    std::string filename = "out.txt";
    int mIter = 300;
    WriteResultsToFile(filename, x, BasicMatrix, ValuesInBasicNodes,
                       tau, N, mIter, dq, q, mem);
//============================================================


    // std::vector<double> Matrix(N - 1);
    // std::vector<double> mem1(N - 1);
    // std::vector<double> B(N - 1);
    // std::vector<double> b(N - 1);


    //double epss = BSolver(x, Matrix, B, b, tau, 4, mIter, mem, mem1);


    return 0;
}