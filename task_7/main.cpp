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


    Fourier(x, p, ValuesInBasicNodes);
    Fullx[0] = Nodes[0];
    for (int i=0;i<N-1;++i){
        Fullx[i+1] = x[i];
    }
    Fullx[N] = Nodes[N];


 
    // std::cout<<"check: "<<std::endl; //<- correct
    // printVector(MultiplyMatrixByVector(FullMatrix, Fullx));
    // printVector(ValuesInNodes);
 
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

    // for (int i = 0; i < N-1; ++i) {
    //     for (int j = 0; j < N-1; ++j) {
    //         if (i==j) MatrixA[i * (N-1) + j] = 1.0; 
    //         else MatrixA[i * (N-1) + j] = 0.0;   
    //     } 
    // }

    // for (int i = 0; i < N-1; ++i) {
    //     b[i] = 1.0;   
    // }
    //printMatrix(MatrixA);
    tau = 1.;


    double pp = BasisMatrixFillWithVariableP(BasicMatrixA); //A
    BasisMatrixFill(pp, BasicMatrixB); //B

    double resid0 = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInNodes, tau, 1, mem, mem1, p); //for x^0
	// std::cout << "Метод 0 запущен" << std::endl;
    std::cout << "resid0 = " << resid0 << std::endl;
    double resid;

    int numberTest = 100;
    double q0 =SearchQ(BasicMatrixA);    //////
    std::cout << "q0 = " << q0 << std::endl;

    std::ofstream fout("output1.txt");
    if (!fout) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    fout << std::fixed << std::setprecision(15);
    for (int iter = 1; iter < numberTest + 1; iter+=4){
        resid = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInNodes, tau, iter, mem, mem1, p); 
        //std::cout << "for iter = " << iter << " started" << std::endl;
        fout << iter << " " << resid << " " << q * resid0 << "\n";
        q *= q0;
    }
    fout.close();

    //printMatrix(BasicMatrix2);

    // BasicMatrix2[0] = 0.020576527053258749307;
    // BasicMatrix2[2] = 0.01434863095211796329;
    // BasicMatrix2[3] = 0.0089143381431437016634;
    // BasicMatrix2[4] = 0.0042504202354287722786;
    // BasicMatrix2[5] = 0.01434863095211796329;
    // BasicMatrix2[6] = 0.044441784361325365809;
    // BasicMatrix2[7] = 0.027610236461134467572;
    // BasicMatrix2[8] = 0.013164758378572473942;
    // BasicMatrix2[9] = 0.0089143381431437016635;
    // BasicMatrix2[10] = 0.027610236461134467572;
    // BasicMatrix2[11] = 0.047833219462407732429;
    // BasicMatrix2[12] = 0.022807221429568085276;
    // BasicMatrix2[13] = 0.0042504202354287722787;
    // BasicMatrix2[14] = 0.013164758378572473942;
    // BasicMatrix2[15] = 0.022807221429568085276;
    // BasicMatrix2[16] = 0.034420678925094271882;

    
    return 0;
}