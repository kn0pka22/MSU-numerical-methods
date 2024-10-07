#include "task_2.hpp"



int main(int argc, char *argv[]){
    int N,k;

    if (argc<3 || argc>4){ printf("Please enter argc=4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3) || (sscanf(argv[2], "%d", &k)!=1) || ((k==0) && (argc==3))){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out N(>2) k(>0)   \n\
        or     enter: ./a.out N(>2) k(=0) filename  \n";
        return -1;
    }

    double* xk = new double[N+1];
    double* yk = new double[N+1];
    double* fmemory = new double[N+1];
    double* U  = new double[N*N];
    double* D  = new double[N*N];
    double* C  = new double[N*N];
    double* phi  = new double[N+1];

    if ((!xk) || (!U) || (!D) || (!C) || (!phi)){
        std::cout<<"Not enough memory!\n";
	 	return -1;
    }
    FillingNodes(xk, N);
    FillingUMatrix(N, U, xk, f);
    FillingDMatrix(N, D, U, phi);
    FillingCMatrix(N, D, C, fmemory, phi);
    

    //for (int i=0;i<N+1;++i) std::cout<<xk[i]<<std::endl;
    // for (int k = 0; k < N * N; k++){
    //     if (k % N == 0) printf("\n");
    //     printf("%20.15lf ", *(D + k));
    // }
    // printf("\n");

    // if (k) WriteToConsole(N, xk, yk, cn, phi);
    // else {
    //     std::string filename = argv[3];
    //     WriteToFile(filename, N, xk, yk, cn, phi); 
    // }


    std::cout<<"In x=2,y=3 f(x,y) = "<<f(2,3)<<"  and Fourier(2,3) = "<<Calc2DFourier(C, N, 2, 3)<<std::endl;

    delete[] xk;
    delete[] yk;
    delete[] fmemory;
    delete[] U;
    delete[] D;
    delete[] C;
    delete[] phi;

	return 0;
}

	