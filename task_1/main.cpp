#include "task_1.hpp"



int main(int argc, char* argv[]){

    int N,k;

    if (argc<3 || argc>4){ printf("Please enter argc=4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3) || (sscanf(argv[2], "%d", &k)!=1) || ((k==0) && (argc==3))){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out N k(>0)   \n\
        or     enter: ./a.out N k(=0) filename  \n";
        return -1;
    }
    //int k = atoi(argv[2]);
    //double h = 1/(N-0.5);

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
    if (k) WriteToConsole(N, xk, yk, cn, phi);
    else {
        std::string filename = argv[3];
        WriteToFile(filename, N, xk, yk, cn, phi); 
    }
   
    delete[] xk; 
    delete[] yk; 
    delete[] phi; 
    delete[] cn; 


    return 0;
}