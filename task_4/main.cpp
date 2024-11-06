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
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
        Please enter:\n\
        N, M, a, b, k\n";
        return -1;
    }

    // Number of nodes required for the algorithm
    // De la Vallée-Poussin polynomial of degree M
	int MM = M + 2;
    


    std::vector<double> Matrix(MM*MM);
    std::vector<double> B(MM);
    //std::vector<double> x(N);   //for what? = Nodes
    std::vector<double> nodes(N);
    std::vector<double> res(MM);
    std::vector<double> ExNodes(3*N-2);
    std::vector<double> values(MM);
    std::vector<double> valuesAll(N);
    std::vector<double> ExValues(3*N-2);
    std::vector<double> ExF(3*N-2);

    std::vector<double> sigma(MM);
    std::vector<int> memory(M+2);


    // std::vector<double> MatrixTest(9);
    // std::vector<double> BTest(3);
    // std::vector<int> memoryTest(9);
    // std::vector<double> resTest(3);
    // MatrixTest[0] = 1;
    // MatrixTest[1] = 1;
    // MatrixTest[2] = 3;
    // MatrixTest[3] = -1;
    // MatrixTest[4] = 1;
    // MatrixTest[5] = 6;
    // MatrixTest[6] = 1;
    // MatrixTest[7] = 1;
    // MatrixTest[8] = 10;
    // BTest[0]=1;
    // BTest[1]=2;
    // BTest[2]=5;
    
    

   

    if (M==N-2){
        //std::cout << "M = N - 2" << std::endl;
        GenerateEquidistantNodes(a, b, nodes);
        for (int i=0;i<MM;++i){
                 sigma[i] = nodes[i]; 
        }
        FillingValues(sigma, values, f, M+2);
        MatrixFill(Matrix,sigma);
        ExtendedNodes(ExNodes, nodes);
        ExtendedValues(ExValues, ExNodes);
        //printVector(ExValues);

        if(solve(Matrix, values, res, memory)){ std::cout<<"smth went wrong, maybe singular matrix\n";}

        std::string filename= "data.txt"; 
        std::ofstream outFile(filename); 
        if (!outFile){
            std::cerr << "Error opening file!" << std::endl;
            return -1;
        }
        outFile <<" | "<< std::setw(10)<< "Nodes"<< std::setw(8)<<" | "<< std::setw(8)<< "Nodes" << std::setw(8)<<"| Погрешность \n";
        outFile << "-------------------------------------------------------\n";

        for (int i = 0; i < MM; i++) {
            if (i == 0){
                outFile << " | " << std::setw(15) << std::setprecision(15) << (i + 1)
                        << " | " << std::setw(15) << std::setprecision(15) << sigma[i]
                        << " | " << std::setw(15) << std::setprecision(15) << res[1] << " \n";
            }
            else {
                outFile << " | " << std::setw(15) << std::setprecision(15) << sigma[i] << " \n";
            }
        }

        outFile.close(); 
        ExtendedF(ExF, res, ExNodes, MM);
        std::cout<<"RESIDUAL: "<<delta(ExValues, ExF, MM)<<std::endl;
    


        std::string filename2 = "graph.txt"; 
        std::ofstream outFile2(filename2); 
        if (!outFile2){
            std::cerr << "Error opening file!" << std::endl;
            return -1;
        }
        for (int i = 0; i < MM - 1; i++) {
            outFile2 << std::setprecision(15) << res[i] << std::endl; 
        }
        outFile2.close(); 



        std::string filename3 = "ans.txt"; 
        std::ofstream outFile3(filename3); 
        outFile3 << " N | X | Истинное значение | Приближенное значение | Разность \n";
        outFile3 << "------------------------------------------------------------------------------------------------------\n";

    
        for (int i = 0; i < 3 * MM - 2; i++) {
            outFile3 << " | " << std::setw(15) << std::setprecision(15) << (i + 1)
                    << " | " << std::setw(15) << std::setprecision(15) << ExNodes[i]
                    << " | " << std::setw(15) << std::setprecision(15) << ExValues[i]
                    << " | " << std::setw(15) << std::setprecision(15) << ExF[i]
                    << " | " << std::setw(15) << std::setprecision(15) << (ExValues[i] - ExF[i]) << " \n";
        }

        outFile3.close(); 
    }
    else if(M<N-2){
        int cnt = 1;
        GenerateEquidistantNodes(a, b, nodes);
        CreateSigma(sigma, nodes, MM, N);
        FillingValues(sigma, values, f, M+2);
        FillingValues(nodes, valuesAll, f, N);
        MatrixFill(Matrix,sigma);
        ExtendedNodes(ExNodes, nodes);
        ExtendedValues(ExValues, ExNodes);
        if(solve(Matrix, values, res, memory)){ std::cout<<"smth went wrong, maybe singular matrix\n";}

        std::string filename = "data.txt"; 
        std::ofstream outFile(filename); 
        if (!outFile){
            std::cerr << "Error opening file!" << std::endl;
            return -1;
        }
        outFile <<" | "<< std::setw(10)<< "Nodes"<< std::setw(8)<<" | "<< std::setw(8)<< "Nodes" << std::setw(8)<<"| Погрешность \n";
        outFile << "-------------------------------------------------------\n";

        for (int i = 0; i < MM; i++) {
            if (i == 0){
                outFile << " | " << std::setw(15) << std::setprecision(15) << (i + 1)
                        << " | " << std::setw(15) << std::setprecision(15) << sigma[i]
                        << " | " << std::setw(15) << std::setprecision(15) << res[1] << " \n";
            }
            else {
                outFile << " | " << std::setw(15) << std::setprecision(15) << sigma[i] << " \n";
            }
        }

        outFile.close(); 
        cnt+=1;
        ExtendedF(ExF, res, ExNodes, MM);

        bool finish = 0;
        while(!finish){
            
        }
        


    }


    
            

 


    return 0;
}
