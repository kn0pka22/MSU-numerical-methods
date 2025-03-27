#include "funcs.hpp"



int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <N> <p>\n"
                  << "N   – number of grid nodes (N >= 1)\n"
                  << "p   – parameter (p > 0)\n";
        return -1;
    }


    int N;
    double p;
    if ((sscanf(argv[1], "%d", &N) != 1) || (N < 1) || (sscanf(argv[2], "%lf", &p) != 1) || (p < 0)) {
        std::cerr << "Invalid input! N must be >= 1 and p must be > 0.\n";
        return -1;
    }

    const int numOfTests = 2;
    std::vector<double> errors(numOfTests);
    std::vector<double> h(numOfTests);

    ClearFile("error_data.txt");

    
    int NCopy = N;
    // N = N/2;

    for (int test = 0; test < numOfTests; ++test) {
        std::vector<double> BasicNodes(N - 1);
        std::vector<double> ValuesInBasicNodes(N - 1);
        std::vector<double> x(N - 1);
        std::vector<double> BasicMatrix((N - 1) * (N - 1));
        std::vector<double> mem(N - 1);

        FillingBasicNodes(BasicNodes);
        FillingValues(BasicNodes, ValuesInBasicNodes, f);
        BasisMatrixFill(p, BasicMatrix);

        Fourier(x, p, ValuesInBasicNodes);

        errors[test] = ErNorm(BasicMatrix, ValuesInBasicNodes, x, N, mem);
        h[test] = 1.0 / (N - 0.5);

        if (test != 0){
            double convergenceRate = log(errors[test - 1] / errors[test]) / log(h[test - 1] / h[test]);
            WriteErrorToFile(N, fabs(convergenceRate), "error_data.txt");
        }

        N *= 2;
    }

    double h1 = 1.0 / (NCopy   - 0.5);
    double h2 = 1.0 / (NCopy*2 - 0.5);
    std::cout<<"N1 = "<<NCopy/2<<std::endl;
    std::cout<<"N2 = "<<NCopy<<std::endl;



    double a = log(errors[0] / errors[1]);
    double b = log(h1 / h2);
    double p_calculated = fabs(a / b);
    // "(h1/h2)^p = err1/err2: p = " 
    std::cout << " p =  "<< p_calculated << std::endl;

 
    return 0;
}