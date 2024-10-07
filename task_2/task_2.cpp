#include "task_2.hpp"


double f(double x, double y){
    return x * (1 - x) * y * (1 - y) * cos(x * x) * cos(y * y);
}


void FillingNodes(double* xk, int N){

    double h = 1/(N-0.5);

    xk[0] = -h/2.;
    for (int i=1;i<N+1;++i){
        xk[i] = xk[i-1]+h;
    }

}

void FillingUMatrix(int N, double* U, double* xk, double (*f)(double, double)){
    for (int i = 1; i < N + 1; ++i){
        for (int j = 1; j < N + 1; ++j){
            if (i == 1 || i == N || j == 1 || j == N)
                U[(i - 1) * (N) + (j - 1)] = 0;
            else
                U[(i - 1) * (N) + (j - 1)] = f(xk[i - 1], xk[j - 1]);
        }
    }
}

double PhiCalculate(int n, int k, double h){
    return sin(M_PI * n * (-h/2. + k*h));
}

void PhiVectorCalculate(int N, int n, double* phi){  //phi_k^{(n)} for fixed k
    double h = 1/(N-0.5);
    //phi[0] = 0;
    for (int k=1;k<N;++k){
        phi[k] = PhiCalculate(n, k, h);
        
    }
    //phi[N]=0;
}

double ScalarProduct(double* ar1,double* ar2, int N){
    double res = 0;
    for (int k=1;k<N; ++k){
        
        res+=ar1[k] * ar2[k];   //the original scalar product requires multiplication by h, but it cancels out
        //std::cout<<res<<std::endl;
    }
    return res;
}

void CoeffCalculate(int N, double* yk, double* phi, double* cn){
    double a = 0;
    double b = 0; 
    for (int n=1; n<N; ++n){
        PhiVectorCalculate(N, n, phi);
        a = ScalarProduct(yk, phi, N);
        b = ScalarProduct(phi,phi, N);
        cn[n] = a/b;
         std::cout<<*(cn + n * N)<<std::endl;
        
    }
}
void FillingDMatrix(int N, double* D, double* U, double* phi){
    for (int j = 0; j < N; ++j){
        D[j * N] = 0;
        D[j * N + N - 1] = 0;
        CoeffCalculate(N - 1, U + j * N, phi, D + j * N);
    }
}

void FillingCMatrix(int N, double* D, double* C, double *fmemory, double* phi){
    for (int i = 0; i < N; ++i){
        C[i * N] = 0;
        C[i * N + N - 1] = 0;
        for (int h = 0; h < N; ++h){
            fmemory[h] = D[h * N + i];
        }
        CoeffCalculate(N - 1, fmemory, phi, C + i * N);
       
    }
}

double Calc2DFourier(double* C, int N, double x, double y){
    double res = 0;
    for (int m = 1; m < N + 1; ++m){
        for (int n = 1; n < N + 1; ++n){
            res += C[(m - 1) * N + (n - 1)] * sin(M_PI * (m - 1) * x) * sin(M_PI * (n - 1) * y);
        }
    }
    return res;
}



// void WriteToConsole(int N, double* xk, double* yk, double* cn, double* phi){
//     double h = 1/(N-0.5);
//     std::cout<< "         xk                   yk                      yk*              "<<std::endl;
//     for (int i = 1; i < (N+1); ++i){
//         yk = 
//         //CoeffCalculate(N, yk, phi, cn);
//         //FourierCompute(cn, N, (- h/2. + i* h));
//         std::cout << std::setprecision(15) << std::fixed \
//         << std::setw(20) << xk[i] << " " \
//         << std::setw(20) << yk[i] << " " \
//         << std::setw(20) << Calc2DFourier(cn, N, xk[i]) << std::endl;
        
//     }
// }

// void WriteToFile(const std::string& filename, int N, double* xk, double* yk, double* cn, double* phi) {
//     double h = 1/(N-0.5);
//     std::ofstream outFile(filename);
//     if (outFile.is_open()) {
//         outFile << "         xk                   yk                      yk*              "<<std::endl;
//         for (int i = 1; i < N+1; ++i){
//             //CoeffCalculate(N, yk, phi, cn);
//             //FourierCompute(cn, N, (- h/2. + i* h));
//             outFile << std::setprecision(15) << std::fixed \
//             << std::setw(20) << xk[i] << " " \
//             << std::setw(20) << yk[i] << " " \
//             << std::setw(20) << Calc2DFourier(cn, N, xk[i]) << std::endl;
//         }
//         outFile.close();
//         std::cout << "Information successfully written to file!" << std::endl;
//     } 
//     else {
//         std::cerr << "Error opening file" << std::endl;
//     }


//}