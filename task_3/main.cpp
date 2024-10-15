#include "task_3.hpp"



int main(int argc, char* argv[]){
    int N,k;

    if (argc<3 || argc>5){ printf("Please enter argc 4 or 5!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<2) || (sscanf(argv[2], "%d", &k)!=1)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter:\n\
        to write to console         enter: ./a.out N(>2) k(>2)   \n\
        to write to file Pn         enter: ./a.out N(>2) k(=0) filename  \n\
        to write to file Ln         enter: ./a.out N(>2) k(=1) filename  \n\
        to write to file Ln and Pn  enter: ./a.out N(>2) k(=2) filename1 filename2  \n";
        return -1;
    }
    double point = 2.2;

    std::vector<double> M(N*N);
    std::vector<double> B(N) ;
    std::vector<int> memory(N);
    std::vector<double> coeffs(N); 
    std::vector<double> nodes(N); 
    std::vector<double> values(N); 

    double a = -1., b= 1.;

    //-------------------calculating P_n----------------------------//
    GenerateEquidistantNodes(a, b, f, nodes, values);
    //GenerateChebyshevNodes(a, b, f, nodes, values);
    MatrixFill(M, nodes);
    VecFill(B, values);
    // std::cout<<"----------------------------Nodes----------------------------"<<std::endl; printVector(nodes);
    
    // std::cout<<"----------------------------Matrix----------------------------"<<std::endl; printMatrix(M);
    
    // std::cout<<"----------------------------Values----------------------------"<<std::endl; printVector(B);

    if (solve(M, B, coeffs, memory)== 0) {
        // std::cout << "---------------------------Solution---------------------------" << std::endl;
        // for (double val : coeffs) {
        //     std::cout << std::setprecision(4) << val << " ";
        // }
        // std::cout << std::endl;
        std::cout<<"Result of Pn calculation:  "<<PnCalculation(coeffs, point)<<std::endl;
        
    } else {
        std::cout << "Singular matrix or error in solving." << std::endl;
    }


    //-------------------calculating L_n----------------------------//
    LnCalculation(nodes, values, point);
    std::cout<<"Result of Ln calculation:  "<<LnCalculation(nodes, values, point)<<std::endl;



    std::string filename = argv[3]; 
    std::string filename2;

    switch (k) {
        case 0:
            WriteToFilePn(a, b, filename, coeffs);
            break;
        case 1:
            WriteToFileLn(a, b, filename, nodes, values);
            break;
        case 2:
            WriteToFilePn(a, b, filename, coeffs);
            filename2 = argv[4];
            WriteToFileLn(a, b, filename2, nodes, values);
            break;
        default:
            std::cerr << "no file to write :( " << std::endl;
            break;
    }




    

    //printVector(b);
    // std::vector<double> M = {
    //     2, 1, -1,  // Row 1
    //     -3, -1, 2, // Row 2
    //     -2, 1, 2   // Row 3
    // };
    // std::vector<double> b = {8, -11, -3}; 
    

    

    return 0;
}

