#include "task_4.hpp"

int main (int argc, char *argv[]){

    //N - num of nodes
    //M - degree of polynomial
    int N,M;
    double a,b;
    int k;
    if (argc<6 || argc>7){ std::cout<<"Please enter argc 6 or 7!\n"; return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (sscanf(argv[2], "%d", &M) != 1) ||
     (N<2) || (sscanf(argv[3], "%lf", &a)!=1) || (sscanf(argv[4], "%lf", &b)!=1) ||
     (sscanf(argv[5], "%d", &k)!=1)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * M - degree of polynomial, \n \
        * [a,b] - segment, \n \
        * k - any value or 0 for file output, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter:\n\
        N, M, a, b, k\n";
        return -1;
    }

    // Number of nodes required for the algorithm
    // De la Vallée-Poussin polynomial of degree M
	int MM = M + 2;
    


    std::vector<double> Matrix(MM*MM);
    std::vector<double> nodes(N);
    std::vector<double> res(MM);
    std::vector<double> ExNodes(3*N-2);
    std::vector<double> values(MM);
    std::vector<double> valuesAll(N);
    std::vector<double> ExValues(3*N-2);
    std::vector<double> ExF(3*N-2);
    std::vector<double> sigma(MM);
    std::vector<int> memory(M+2);

    
    std::string filename;
    if (M==N-2){
        
        //GenerateEquidistantNodes(a, b, nodes);
        GenerateChebyshevNodes(a,b,nodes);
        for (int i=0;i<MM;++i){
                sigma[i] = nodes[i]; 
        }
        FillingValues(sigma, values, f, M+2);
        MatrixFill(Matrix,sigma);
        
        ExtendedNodes(ExNodes, nodes);
        ExtendedValues(ExValues, ExNodes);
        //printVector(ExValues);

        if(solve(Matrix, values, res, memory)){ std::cout<<"smth went wrong, maybe singular matrix\n";}
        //printVector(res);
        

        filename= "data.txt"; 
        WriteToFile(a, b, filename, res);

        ExtendedF(ExF, res, ExNodes, MM);
        std::cout<<"RESIDUAL: "<<delta(ExValues, ExF, MM)<<std::endl;
    
    }
    else if(M<N-2){
        GenerateEquidistantNodes(a, b, nodes);
        //printVector(nodes);
        CreateSigma(sigma, nodes, MM, N);
        //printVector(sigma);
        FillingValues(sigma, values, f, MM);
        FillingValues(nodes, valuesAll, f, N);
        MatrixFill(Matrix,sigma);
        //printMatrix(Matrix);
        //print();


        if(solve(Matrix, values, res, memory)){ std::cout<<"smth went wrong, maybe singular matrix\n";}
        FillingValues(sigma, values, f, MM);
       

        bool flag = 1;
        int  iter = 1;
        while((flag) && (iter < MAX_ITERATIONS)){
           //std::cout<<"HERE! and num of iteration = "<<iter<<std::endl;
            flag = MaxDeviation(nodes, sigma, res, values, valuesAll, MM, N);
            //flag=0;
            if(flag){
                MatrixFill(Matrix,sigma);
                if(solve(Matrix, values, res, memory)){ std::cout<<"smth went wrong, maybe singular matrix\n";}
                FillingValues(sigma, values, f, MM);
                iter += 1;
            }
        }
        //filename= "data.txt"; 
        //WriteToFile(a, b, filename, res);

        ExtendedNodes(ExNodes, nodes);
        ExtendedValues(ExValues, ExNodes);
        ExtendedF(ExF, res, ExNodes, MM);
        std::cout<<"RESIDUAL: "<<delta(ExValues, ExF, MM)<<std::endl;

   
        
    }

    return 0;
}









