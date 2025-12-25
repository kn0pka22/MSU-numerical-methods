#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "funcs.hpp"

const int test_for_conv = 10;

int main() {
    int N = 1000;
    double A = 1000.;
    double y0 = 1.;
    double y1, y2, y3, y4;

    std::ofstream outFile("out.txt");
    if(outFile.is_open()) {
        for(int i = 1; i <= test_for_conv; i++) {
            y1 = scheme_1(y0, A, N);
            y2 = scheme_2(y0, A, N);
            y3 = scheme_3(y0, A, N);
            y4 = scheme_4(y0, A, N);
            // outFile << std::setprecision(15) << N << " " 
            //          << log(y1) << " " << log(y2) << " " 
            //          << log(y3) << " " << log(y4) << std::endl;
            outFile << std::setprecision(15) << N << " " 
                     << y1 << " " << y2 << " " 
                     << y3 << " " << y4 << std::endl;
            N *= 2;
        }
    } else { 
        std::cerr << "ERROR: Unable to create out.txt" << std::endl;
        return 1;
    }
    outFile.close();
    
    return 0;
}