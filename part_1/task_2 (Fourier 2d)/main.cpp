#include "task_2.hpp"



int main(int argc, char *argv[]){
    int N,k;

    if (argc<3 || argc>4){ printf("Please enter argc 3 or 4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3) || (sscanf(argv[2], "%d", &k)!=1) || ((k==0) && (argc==3))){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter:\n\
        to write to console enter: ./a.out N(>2) k(>0)   \n\
        to write to file    enter: ./a.out N(>2) k(=0) filename  \n\
        for pcalculation    enter: ./a.out N(>2) k(<0) \n";
        return -1;
    }

    double* xk      = new double[N+1];
    double* fmemory = new double[N+1];
    double* U       = new double[(N+1)*(N+1)];
    double* D       = new double[(N+1)*(N+1)];
    double* C       = new double[(N+1)*(N+1)];
    double* phi     = new double[N+1];
 

    if ((!xk) || (!U) || (!D) || (!C) || (!phi) || (!fmemory)){
        std::cout<<"Not enough memory!\n";
	 	return -1;
    }


    if (k>0){ 
        //std::cout<<"DURATION (fToC)= "<<WriteToConsole(N, xk, U, C, D, fmemory, phi)<<std::endl;
        double err1 = normFunction(f, C, N);
        std::cout << "err1 = " << fabs(err1) << std::endl;
    }
    else if (!k) {
        std::string filename = argv[3];
        WriteToFile(filename, N, xk, U, C, D, fmemory, phi); 
    }
    else pcalculate(N);
    
    
    
    delete[] xk;
    delete[] fmemory;
    delete[] U;
    delete[] D;
    delete[] C;
    delete[] phi;

	return 0;
}

	