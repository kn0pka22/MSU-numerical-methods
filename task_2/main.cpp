#include "task_1.hpp"


void PrintTable(double *nodes, int N, double *coefs, double (*f)(double)){
    std::cout<<"   Nodes               Fourier              f\n";
    for (int i = 1; i < N - 1; ++i)
    {
        double xi = nodes[i];
        double delta = nodes[i+1] - xi;
        delta /= 3.;
        std::cout << std::setprecision(15) << std::fixed \
        << std::setw(20) << xi << " " \
        << std::setw(20) << FourierCompute(coefs, N, xi) << " " \
        << std::setw(20) << f(xi) << std::endl;
        double xi1 = xi+ delta;
        std::cout << std::setprecision(15) << std::fixed \
        << std::setw(20) << xi1 << " " \
        << std::setw(20) << FourierCompute(coefs, N, xi1) << " " \
        << std::setw(20) << f(xi1) << std::endl;
        double xi2 = xi + 2*delta;
        std::cout << std::setprecision(15) << std::fixed \
        << std::setw(20) << xi2 << " " \
        << std::setw(20) << FourierCompute(coefs, N, xi2) << " " \
        << std::setw(20) << f(xi2) << std::endl;
    }
    double xi = nodes[N - 1];
    std::cout << std::setprecision(15) << std::fixed \
        << std::setw(20) << xi << " " \
        << std::setw(20) << FourierCompute(coefs, N, xi) << " " \
        << std::setw(20) << f(xi) << std::endl;
}

void WriteResult(double *nodes, int N, double *coefs, double (*f)(double), const std::string& filename){
    double xi = 0;
    double delta = 0; 
    std::ofstream outfile(filename); 
    if (outfile.is_open()) {
        for (int i = 1; i < N - 1; ++i){
            xi = nodes[i];
            delta = nodes[i+1] - xi;
            delta /= 3.;
            outfile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << FourierCompute(coefs, N, xi) << " " \
            << std::setw(20) << f(xi) << std::endl;
            double xi1 = xi+ delta;
            outfile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi1 << " " \
            << std::setw(20) << FourierCompute(coefs, N, xi1) << " " \
            << std::setw(20) << f(xi1) << std::endl;
            double xi2 = xi + 2*delta;
            outfile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi2 << " " \
            << std::setw(20) << FourierCompute(coefs, N, xi2) << " " \
            << std::setw(20) << f(xi2) << std::endl;
        }
        outfile.close(); 
        std::cout << "Congratulations! Result in file, you can check them right now!" << filename << std::endl;
    } 
    else{
        std::cerr << "smth went wrong in opening file :( " << filename << std::endl;
    }
}

void WriteSkrypt(int N, const std::string& testfilename, const std::string& filename){
    
    std::ofstream OutFile(filename); 
    std::ofstream OutTestFile(testfilename); 
    if (OutFile.is_open() && OutTestFile.is_open()) {

        OutFile<< "#! /usr/bin/gnuplot -persist\n";
        OutFile<< "set terminal png size 1000,1000 enhanced font \"Helvetica Bold, 20\"\n";
        OutFile<< "set output \""<<testfilename<<".png\"\n\n";

        OutFile<< "set style line 1 lt 1 linecolor rgb \"red\" lw 1 pt 1\n";
        OutFile<< "set style line 2 lt 1 linecolor rgb \"blue\" lw 1 pt 1\n";

        //    fprintf(out, "set yrange [0:5]");
        OutFile<< "set xrange [0:1]\n";

        OutFile << "set title \"" << testfilename << " - " << N << " nodes \"\n";

        OutFile << "set grid\n\n";

        OutFile<< "plot  \"" << testfilename << "\" using 1:2 ls 1 title \"Interpolation Fourier Row\" with lines, "; \
        OutFile << "\"" << testfilename << "\" using 1:3 ls 2 title \"Original function\", ";
        OutFile.close(); OutTestFile.close(); 
        std::cout << "Congratulations! Result in file, you can check them right now!" << filename << std::endl;
  
    }
    else{
        std::cerr << "smth went wrong in opening file :( " << filename << std::endl;
    }

}



int main(int argc, char *argv[]){

	int N; //количество узлов сетки
	//
    

    double* x;
    double* y;
    double* phi;  //базис
    double* c;    //искомые коэффициенты

	if (argc<3 || argc>4){ printf("Please enter argc=4!\n"); return -1;} 
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3)){ 
        std::cout<<"Invalid input!\n \
        * N – number of grid nodes, \n \
        * filename – the name of the file to which the result should be written. This argument is missing if k! = 0.\n\n\
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

    N=atoi(argv[1]);

    if(CoeffCalculate(N, c, f, x, y, phi)) {
        std::cout<<"c_m not calculated - error"<<std::endl;
        return -1;
    }
    std::cout<<"Great! Fourier series coefficients calculated!"<<std::endl;
    
    //std::ofstream f_out("out.txt"); 
    int k=atoi(argv[2]);
    if (k == 0){
        std::cout<<"no files to write answer"<<std::endl;
        //return 0;
    }
    else{ std::cout<<"file opened!"<<std::endl;}

    if (k!=0) PrintTable(x, N, c, f);
    else{
        //filename = 
        //std::ofstream outfile("printAll.gpi"); 
        //if (outfile.is_open()) {
            std::string filename = argv[3];
            std::string filename2 = "printAll.gpi";
            WriteSkrypt(N, filename, filename2);
            //WriteResult(x, N, c, f, filename);
        //}
        //outfile.close(); 
    }

    std::cout<<"Congratulations! Result in file, you can check them right now!"<<std::endl;

	free(x);
    free(y);
    free(phi);
    free(c); 

	return 0;
}

	