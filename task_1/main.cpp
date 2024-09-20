#include "task_1.hpp"



int main(int argc, char* argv[]){

    int N;

    if (argc<3 || argc>4){ printf("Please enter argc=4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out N k(>0)   \n\
        or     enter: ./a.out N k(=0) filename  \n";
        return -1;
    }
    double h = 1/(N-0.5);

    double* xk  = new double[N+1];
    double* yk  = new double[N+1];
    double* phi = new double[N+1];
    double* cn  = new double[N+1];


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
    

    std::cout<<"      xk         yk         yk*              "<<std::endl;
    for (int i = 1; i < N; ++i){
        CoeffCalculate(N, i, yk, phi, cn);
        FourierCompute(cn, N, (- h/2. + i* h));
        std::cout << std::setprecision(5) << std::fixed \
        << std::setw(10) << xk[i] << " " \
        << std::setw(10) << yk[i] << " " \
        << std::setw(10) << FourierCompute(cn, N, xk[i]) << std::endl;
        
    }

    delete[] xk; 
    delete[] yk; 
    delete[] phi; 
    delete[] cn; 


    return 0;
}