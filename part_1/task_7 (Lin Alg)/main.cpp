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

    std::vector<double> MatrixA((N - 1) * (N - 1));
    std::vector<double> b(N - 1);

    

    std::vector<double> BasicMatrix((N - 1) * (N - 1));
    std::vector<double> BasicMatrixA((N - 1) * (N - 1));
    std::vector<double> BasicMatrixB((N - 1) * (N - 1));
    std::vector<double> FullMatrix((N + 1) * (N + 1));

    std::vector<double> x(N - 1);
    std::vector<double> Fullx(N + 1);

    std::vector<double> mem(N - 1);
    std::vector<double> mem1(N - 1);



    FillingNodes(Nodes);
    FillingBasicNodes(BasicNodes);
    FillingValues(Nodes, ValuesInNodes, f);
    FillingValues(BasicNodes, ValuesInBasicNodes, f);

 
    BasisMatrixFill(p, BasicMatrix);
    FullMatrixFill(p, FullMatrix);
    //printMatrix(FullMatrix);

//============================TASK 1================================

    Fourier(x, p, ValuesInBasicNodes);
    Fullx[0] = Nodes[0];
    for (int i=0;i<N-1;++i){
        Fullx[i+1] = x[i];
    }
    Fullx[N] = Nodes[N];


 
    // std::cout<<"check: "<<std::endl; //<- correct
    // printVector(MultiplyMatrixByVector(FullMatrix, Fullx));
    // printVector(ValuesInNodes);


    std::cout<<"check: "<<std::endl; //<- correct
    printVector(MultiplyMatrixByVector(BasicMatrix, x));
    printVector(ValuesInBasicNodes);
 
//============================TASK 2================================
    
    double EigenValueMin = Lambda(1, N, p);
    double EigenValueMax = Lambda(N-1, N, p);
    
    double tau = 2. / (EigenValueMin + EigenValueMax);
    double q = (EigenValueMax - EigenValueMin) / (EigenValueMax + EigenValueMin);
    double dq = 1;


    std::string filename = "out.txt";
    int mIter = 300; 
    WriteResultsToFile(filename, x, BasicMatrix, ValuesInBasicNodes,
                       tau, N, mIter, dq, q, mem);
//============================TASK 3================================


    std::ofstream fout("output1.txt");
    if (!fout) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    double pp = BasisMatrixFillWithVariableP(BasicMatrixA); //A
    BasisMatrixFill(pp, BasicMatrixB); //B

    //std::cout<<"pp = "<<pp<<std::endl;
     int numberTest = 300;
    //double q0 =SearchQ(BasicMatrixA);
    double q0 = q;
    double resid0 = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInBasicNodes, tau, 1, mem, mem1, pp); //for x^0
    //std::cout << "resid0 = " << resid0 << std::endl;
    fout << "0  " << resid0 << " " << q0 * resid0 << "\n";
    double resid;

    //printMatrix(BasicMatrixA);
    //printMatrix(BasicMatrixB);
    
   
    //q = q0;    
    //std::cout << "q0 = " << q0 << std::endl;

    // for(int i = 0; i < N-1; i++)
    // {
    //     ValuesInBasicNodes[i] = 0.;
    // }
    // ValuesInBasicNodes[N/2] = 1.;
    
    tau = 1.;
    fout << std::fixed << std::setprecision(15);
    for (int iter = 0; iter < numberTest+2; iter+=1){
        resid = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInBasicNodes, tau, iter, mem, mem1, pp); 
        //std::cout << "for iter = " << iter << " started" << std::endl;
        fout << iter << " " << resid << " " << q0 * resid0 << "\n";
        q0 *= q0;
    }
    fout.close();
    
    return 0;
}
