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

    std::vector<double> A((N + 1) * (N + 1), 0.0);
    std::vector<double> b(N + 1, 0.0);
    std::vector<double> x(N + 1, 1.);
    std::vector<double> bCheck(N + 1, 0.0);
    std::vector<double> mem(N + 1, 0.0);

   
    MatrixFill(p, A);
    //printMatrix(A);
    RightSideFill(A, b);
    //b[0]=b[N]=0; 
    //printVector(b);

    FourierMethod(x, N, p, b);
    
    // //check:
    // bCheck = MultiplyMatrixByVector(A, x);
    // printVector(b);
    // printVector(bCheck);  //correct

    
    double m = Lambdan(1, N + 1, p);
    double M = Lambdan(N, N + 1, p);
    double tau = 2. / (m + M);

    double q = (M - m) / (M + m);
    double dq = 1;


    std::string filename = "inp.txt";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading!" << std::endl;
        return 1;
    }
    file<< std::left 
        << std::setw(15) << "N" 
        << std::setw(20) << "Eps" 
        << std::setw(20) << "dq*fnorm" 
        << std::endl;
    double fnorm = Richardson(x, A, b, tau, N, 0, mem);
    double eps = 0.;
    for (int n = 0; n < m; n += 1){
        eps = Richardson(x, A, b, tau, N, n, mem);
        file << std::left 
             << std::setw(15) << std::fixed << std::setprecision(10) << n
             << std::setw(20) << std::fixed << std::setprecision(10) << eps
             << std::setw(20) << std::fixed << std::setprecision(10) << dq * fnorm
             << std::endl;
        dq *= q;
    }
    return 0;
}