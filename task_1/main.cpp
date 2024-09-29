#include "task_1.hpp"



int main(int argc, char* argv[]){

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

    //int Numknots = 2;

    double* xk        = new double[N+1];
    double* yk        = new double[N+1];
    double* phi       = new double[N+1];
    double* cn        = new double[N+1];

    // double* cnForP    = new double[N+1];
    // double* cnForP2   = new double[2*N+1];
    // double* xkForP2   = new double[2*N+1];
    // double* ykForP2   = new double[2*N+1];
    // double* phiForP  = new double[2*N+1];
    // double* phiForP2  = new double[2*N+1];
    
    //double* lognorms  = new double[Numknots];
    //double* loghs     = new double[Numknots];


    if ((!xk) || (!yk) || (!phi) || (!cn)){
        std::cout<<"Not enough memory!\n";
	 	return -1;
    }

    //xk = nullptr;    //    
    FillingNodes(xk, N);
    try {
        FillingValues(xk, yk, f, N);
    } 
    catch (const std::exception& e){
        std::cerr << "Error: " << e.what() << std::endl;
    }
    if (k) WriteToConsole(N, xk, yk, cn, phi);
    else {
        std::string filename = argv[3];
        WriteToFile(filename, N, xk, yk, cn, phi); 
    }



    // for (int i = 1; i < N+1; ++i){
    //     //std::cout<<cnForP[i]<<std::endl;
    //     CoeffCalculate(N, i, yk, phiForP, cnForP);
    //     std::cout<<cnForP[i]<<std::endl;
    // }
    // // for (int i = 1; i < 2*N+2; ++i){
    //     CoeffCalculate(2*N, i, yk, phi, cnForP2);
    // }

    // double err1 = NormFunction(f,cnForP,  N);
    // double err2 = NormFunction(f,cnForP2,2*N);

    // double h1 = 1/(double)N;
    // double h2 = 1/(double)(2*N);

    // double a = log(err1/err2);
    // double b = log(h1/h2);


    //std::cout<<"(h1/h2)^p = err1/err2:"<<"  ";
    //std::cout<<"p =  "<<a/b<<std::endl;


    // PCalculate(N, cn, lognorms, loghs, Numknots);
    // double p=0;
    // for (int i=0;i<Numknots;i+=2){
    //     p+=fabs((lognorms[i]-lognorms[i+1])/(loghs[i]-loghs[i+1]));
    // }
    // p *= 1/((double)Numknots);
    // std::cout<<"the value of p: "<<p<<std::endl;

    pcalculate(N);

   
    delete[] xk; 
    delete[] yk; 
    delete[] phi; 
    delete[] cn; 

    // delete[] cnForP; 
    // delete[] cnForP2; 
    // delete[] xkForP2; 
    // delete[] ykForP2; 
    // delete[] phiForP; 
    // delete[] phiForP2; 
    //delete[] lognorms; 
    //delete[] loghs; 


    return 0;
}

