#include "funcs.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <fstream>
 

double func(double t, double x) { 
    return exp(-M_PI*M_PI*t)*(-M_PI * M_PI *sin(x)*(x-1) + sin(x)*(x-1)-2*cos(x));  
    // return exp(-M_PI*M_PI*t)*(-M_PI*M_PI*(x*(x-1.))-2);
}
 
double u0(double x) { 
    // return x*(x-1);
    return sin(x)*(x-1);

}

void printVec(double* vec, int N){
    for (int i=0;i<N;++i){
        std::cout << std::fixed << std::setprecision(5) << vec[i] << " ";
    }
    std::cout<<std::endl;
}

double uFunc(double t, double x){
    // return exp(-M_PI*M_PI*t)*x*(x-1);
    return exp(-M_PI*M_PI*t) * sin(x) * (x - 1);


}

void printMatr(double** matr, int N, int M){
    for (int i=0;i<N;++i){
        for (int j=0;j<M;++j){
            std::cout << std::setw(10) << std::fixed << std::setprecision(2) << matr[i][j];
        }
        std::cout<<std::endl;
    }
    // std::cout<<std::endl;
}


int findSteps(double T, double tau) {
    if (tau < 1e-10) {
        throw std::invalid_argument("tau is too small or zero");
    }
    return static_cast<int>(std::floor((T + 1e-10) / tau));  
}



// неявная схема 
void solveImplicit(int M, double levels, double tau, double** u) {
    double h   = 1.0 / (M - 0.5);
    
    u[0][0] = u0(-h/2.);

    // std::cout<<"u[0][0] = "<<u[0][0]<<std::endl;
    for (int j = 1; j < M+1; ++j) {
        u[0][j] = u0(-h/2. + j*h);
        // std::cout<< " -h/2. + j*h = "<<-h/2. + j*h<<std::endl;
        // std::cout<< " sin(M_PI * ("<<-h/2. + j*h<<")) = "<<u0(-h/2. + j*h)<<std::endl;
    }

    // std::cout<<"u[0][M] = "<< u[0][M]<< " and must be = 0 ! "<<std::endl; 

    double* f = new double[M+1];

    // int Nt      = findSteps(T,tau);
    // std::cout<<" num of steps by axis t = "<<Nt<<std::endl;
    // double lastTau = T - tau * Nt; // 
    double t;
    double* upDiag   = new double[(M+1)];
    double* midDiag  = new double[(M+1)];
    double* downDiag = new double[(M+1)];


    createThreeDiag(upDiag, midDiag, downDiag, M, h, tau);
    // printVec(upDiag, (M+1));
    // printVec(midDiag, (M+1));
    // printVec(downDiag, (M+1));


    for (int i = 1; i < levels+1; ++i) { //для каждого уровня решаем задачу
        t = i * tau;
        f[0] = 0.;  
        f[M] = 0.;
        for (int j = 1; j < M; ++j) {
            f[j] = u[i-1][j]/tau + func(t, -h/2. + j*h);
        }
        SolveTridiagonal(upDiag, midDiag, downDiag, f, M);
        for (int j=0;j<M+1;++j){
            u[i][j] = f[j];
        }

    }
    // if (lastTau>1e-10){
    //     std::cout<<" lastTau>1e-10 \n ";
    //     f[0] = 0.;  //exactly 0?
    //     f[M] = 0.;
    //     for (int j = 1; j < M; ++j) {
    //         f[j] = u[Nt][j]/tau + func(t, -h/2. + j*h);
    //     }
    //     createMatrix(threeDiag, M, h, tau);
    //     SolveTridiagonal(threeDiag, f, M);
    //     for (int j=0;j<M+1;++j){
    //         u[Nt+1][j] = f[j];
    //     }
    // }

    delete[] f;
    delete[] upDiag;
    delete[] midDiag;
    delete[] downDiag;

    
}


void createThreeDiag(double* upDiag, double* midDiag, double* downDiag, int M, double h, double tau){

    upDiag[0]  = 1.;
    midDiag[0] = 1.;
    downDiag[0]= 0.;

    // upDiag[1]   = -1./(h*h);
    // midDiag[1]  =  3./(h*h) + 1./tau;
    // downDiag[1] = -1./(h*h);

    for (int i=1;i<M;++i){
        upDiag[i]   = -1./(h*h);
        midDiag[i]  =  2./(h*h) + 1./tau;
        downDiag[i] = -1./(h*h);
    }

    upDiag[M]  = 0.;
    midDiag[M] = 1.;
    downDiag[M]= 0.;
}


void SolveTridiagonal(double* upDiag, double* midDiag, double* downDiag, double* rhs, int M){                   

    double* alpha = new double[M + 1];
    double* beta  = new double[M + 1];
    

    alpha[0] = -upDiag[0] / midDiag[0];
    beta[0] = rhs[0] / midDiag[0];

    for (int i = 1; i < M; ++i) {
        double denominator = midDiag[i] + downDiag[i] * alpha[i - 1];
        alpha[i] = -upDiag[i] / denominator;
        beta[i] = (rhs[i] - downDiag[i] * beta[i - 1]) / denominator;
    }

    rhs[M] = (rhs[M] - downDiag[M] * beta[M - 1]) 
             / (midDiag[M] + downDiag[M] * alpha[M - 1]);

    for (int i = M - 1; i >= 0; --i) {
        rhs[i] = alpha[i] * rhs[i + 1] + beta[i];
    }

    delete[] alpha;
    delete[] beta;
}


double error(double** u, int M, int N, double h, double T){
    double sum =0.0;
    for (int i=0;i<M+1;++i){
        sum+=(u[N][i]-uFunc(T,-h/2. + i*h))*(u[N][i]-uFunc(T,-h/2. + i*h));
    }
    sum = sum*h;

    return sqrt(sum);
}

//явная схема
void solveExplicit(int N, int M, double tau, double h, double** u) {
   
    for (int j = 0; j < M+1; ++j) {
        u[0][j] = u0(-h/2. + j*h);
    }
    // std::cout<<"u_m^0 = ";
    // printVec(u[0], M+1);

    // for (int i = 1; i < N+1; ++i) {
    //     // u[i][0] = 0.0;
    //     u[i][M] = 0.0;
    // }

    double r = tau / (h * h);

    double current_t;

    // 1/N <= 1/(2.*(M-0.5)^2)
    
    // Проверка устойчивости
    // if (r >= 0.5) {
    if (N < 2.*(M-0.5)*(M-0.5)){
        std::cerr << "Warning: Explicit scheme may be unstable (tau / (h * h) = " << r << " >= 0.5)\n";
    }

    for (int i = 0; i < N; ++i) {
        // Правое граничное условие: u(n, M) = 0
        u[i+1][M] = 0.0;

        current_t = i * tau;
        for (int j = 1; j < M; ++j) {
            u[i+1][j] = u[i][j] + r * (u[i][j-1] - 2.0 * u[i][j] + u[i][j+1]) 
                       + tau * func(current_t, -h/2. + j*h);
        }
        // u[i+1][0] = tau * (2/(h*h)*(u[i][1]+u[i][0])  + 0*uFunc(current_t,0)) + u[i][0];
        // Левое граничное условие: (u(n,0) + u(n, h/2)) / 2 = δ(n)
        // u[i+1][0] = u[i][0] + tau * ( (2.0/(h*h)) * (u[i][1] + u[i][0]) 
                    // + 2.0 * func(current_t, 0.0) );
        // u[i+1][0] = tau * ((u[i][0]+u[i][1])/(h*h) + uFunc(current_t,0)) + u[i][0];
        u[i+1][0] = - u[i+1][1];
    }
    // std::cout<<"u_m^1 = ";
    // printVec(u[1], M+1);

}

void freeSolution(double** u, int N) {
    for (int i = 0; i < N+1; ++i) {
        delete[] u[i];
    }
    delete[] u;
}

void saveMatrixToFile(double** matr, int rows, int cols, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            outFile << std::setprecision(6) << matr[i][j];
            if (j < cols - 1) outFile << " ";  
        }
        outFile << "\n";  
    }
    outFile.close();
}

void saveVectorToFile(double* vec, int size, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    for (int i = 0; i < size; ++i) {
        outFile << std::setprecision(6) << vec[i] << "\n";
    }
    outFile.close();
}
void runConvergenceTest(const std::string& filename, const std::string& SchemeName, double T, int testCase) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw std::runtime_error("Cannot open file for convergence test results");
    }

    outfile << std::fixed 
            << std::setw(8) << "M"      
            << std::setw(8) << "N"              
            << std::setw(15) << "h"       
            << std::setw(15) << "tau"      
            << std::setw(15) << "error"     
            << std::setw(15) << "log_h"    
            << std::setw(15) << "log_error"   
            << std::endl;
    
    const int test_count = 5;
    const int M_values[] = {10, 20, 40, 80, 160};
    int N_values[test_count];

    // Выбираем N_values в зависимости от testCase
    switch(testCase) {
        case 1:  //  impl M=N
            for (int i = 0; i < test_count; ++i) {
                N_values[i] = M_values[i];
            }
            break;

        case 2:  //  impl M,N -> 2M,4N 
            for (int i = 0; i < test_count; ++i) {
                N_values[i] = 2 * (M_values[i]) * (M_values[i]);
            }
            break;
            
        case 3:  // expl M,N -> 2M,2N 
        {
            int N0 = 2 * (M_values[4]) * (M_values[4]);
            for (int i = 0; i < test_count; ++i) {
                N_values[i] = N0 * (1 << i);  // Умножение на 2^i
            }
            break;
        }
        
        default:
            throw std::invalid_argument("Unknown test case");
    }

    for (int i = 0; i < test_count; ++i) {
        int M = M_values[i];
        int N = N_values[i];
        double h = 1.0 / (M - 0.5);
        double tau = T / N;
        
        double** u = new double*[N+1];
        for (int j = 0; j < N+1; ++j) {
            u[j] = new double[M+1]();
        }

        if (SchemeName == "explicit") {
            solveExplicit(N, M, tau, h, u);
        } else {
            solveImplicit(M, N, tau, u);
        }
        
        double err = error(u, M, N, h, T);
        
        outfile << std::scientific 
                << std::setfill(' ')
                << std::setw(8) << M     
                << std::setw(8) << N       
                << std::setw(15) << h     
                << std::setw(15) << tau   
                << std::setw(15) << err    
                << std::setw(15) << log(h)  
                << std::setw(15) << log(err)
                << std::endl;
        
        freeSolution(u, N);
    }
    
    outfile.close();
    std::cout << "Convergence results saved to " << filename << "\n";
}