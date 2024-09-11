#include "task_1.hpp"

int main(int argc, char *argv[]){

	int N; //количество узлов сетки
	//FILE* in;
    //FILE* out;

    double* x;
    double* y;
    double* phi;  //базис
    double* c;    //искомые коэффициенты

	if (argc<3 || argc>4){ printf("Please enter argc=4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<0)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file from which the matrix should be read. This argument is missing if k! = 0.\n\n\
        Please enter: ./a.out N k(>0)   \n\
        or     enter: ./a.out N k(=0) filename  \n";

        return -1;
    }
	x   = (double*)malloc((N + 1) * sizeof(double));
    y   = (double*)malloc((N + 1) * sizeof(double));
    phi = (double*)malloc((N + 1) * sizeof(double));
    c   = (double*)malloc((N + 1) * sizeof(double));

	if ((!x) || (!y) || (!phi) || (!c)) {
		std::cout<<"Not enough memory!\n";
		return -1;
	}

	free(x);
    free(y);
    free(phi);
    free(c); 

	return 0;
}

	