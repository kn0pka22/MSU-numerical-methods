#include "funcs.hpp"


int NumOfTests = 10;

int main() {
    int N = 1000;
    double A = 1000.;
    double y0 = 1.;
    double y1, y2, y3, y4;

    std::ofstream outFile("data.txt");
    if(!outFile.is_open()) {
        std::cerr << "Error: Could not open data.txt" << std::endl;
        return 1;
    }

    outFile << std::scientific << std::setprecision(15);
    for(int i = 1; i < NumOfTests+1; i++) {
        y1 = scheme_1(y0, A, N);
        y2 = scheme_2(y0, A, N);
        y3 = scheme_3(y0, A, N);
        y4 = scheme_4(y0, A, N);
        
        outFile << log(N) << " " 
                << log(y1) << " " 
                << log(y2) << " " 
                << log(y3) << " " 
                << log(y4) << "\n";
        N *= 2;
    }
    outFile.close();

    std::cout << "Data successfully written to data.txt" << std::endl;
    return 0;
}