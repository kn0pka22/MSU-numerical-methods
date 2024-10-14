#include "task_3.hpp"



int main(int argc, char* argv[]){
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


    std::vector<double> M = {
        2, 1, -1,  // Row 1
        -3, -1, 2, // Row 2
        -2, 1, 2   // Row 3
    };
    std::vector<double> b = {8, -11, -3}; // Right-hand side vector
    std::vector<double> x(N); // Vector to store results
    std::vector<int> memory(N);

    int result = solve(M, b, x, memory);

    if (result == 0) {
        std::cout << "Solution:" << std::endl;
        for (double val : x) {
            std::cout << std::setprecision(4) << val << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Singular matrix or error in solving." << std::endl;
    }

    return 0;
}

