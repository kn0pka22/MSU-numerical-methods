#include "task_7.hpp"

double f(double x) {
    return sin(M_PI * x);
}

void printMatrix(const std::vector<double>& matrix) {
    int n = sqrt(matrix.size());
    for (int i = 0; i < n; ++i) { 
        for (int j = 0; j < n; ++j) { 
            std::cout << std::setw(10) << std::setprecision(4) << matrix[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vec) {
    for (double val : vec) {
        std::cout << std::setw(10) << std::setprecision(4) << val << " ";
    }
    std::cout << std::endl;
}

std::vector<double> MultiplyMatrixByVector(const std::vector<double>& matrix, const std::vector<double>& vec) {
    int n = sqrt(matrix.size());

    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i*n+j] * vec[j];
        }
    }
    return result;
}


void MatrixFill(double p, std::vector<double>& M){
    
    int N = sqrt(M.size());  //N+1
    N--;

    for (int i = 0; i < N+1; i++) {
        for (int j = 0; j < N+1; j++) {
            if (i == j && i > 0 && j > 0 && i < N && j < N){
                M[i * (N + 1) + j] = p + 2.0 * (double)(N * N);
            } 
            else if ((i - j == 1 || i - j == -1) && i > 0 && j > 0 && i < N && j < N){
                M[i * (N + 1) + j] = -(double)(N * N);
            } 
            else{
                M[i * (N + 1) + j] = 0.0;
            }
        }
    }
}

void RightSideFill(const std::vector<double>& A, std::vector<double>& b) {
    
    int N = b.size(); //N+1
    N--; 

    double h = 1./(double)(N-2);
    double xi = h/2.;

    b[0] = 0;
    for (int i = 1; i < N; ++i){
        b[i] = f(xi);
        xi=-h/2. + i*h;
    }
    b[N] = 0;

}



// double psi(int k, int n, int N){
//     return sin(M_M_PI * n * k / (double)(N));
// }


double psi(int k, int n, int N){
    return sin(M_PI * n * k / (double)(N));
}

double Lambdan(int n, int N, double p){
    double lam = p - 2 * N * N * (cos(M_PI * n / (double)(N)) - 1);
    return lam;
}

double Dn(int n, std::vector<double>& f, double p, int N){
    double sp = 0;
    for (int i = 0; i < N; ++i){
        sp += 2 * f[i] * psi(i, n, N) / (double)(N);
    }
    return sp;
}

double FourierMethod(std::vector<double>& y, int N, double p, std::vector<double>& f){
    for (int k = 0; k < N + 1; ++k){
        y[k] = 0;
        for (int n = 1; n < N; ++n){
            y[k] += Dn(n, f, p, N) / Lambdan(n, N, p) * psi(k, n, N);
        }
    }
    return 0.;
}

double Richardson(std::vector<double>& x, const std::vector<double>& A, const std::vector<double>& b, double tau, int n, int mIter, std::vector<double>& mem) {
    for (int k = 0; k < n; ++k){
        x[k] = 0;
        mem[k] = 0;
    }
    for (int m = 0; m < mIter; m++){
        mem = MultiplyMatrixByVector(A, x);
        for (int i = 0; i < n; i++){
            x[i] = x[i] - tau * mem[i] + tau * b[i];
        }
    }
    return ErNorm(A, b, x, n, mem);
}

double ErNorm(const std::vector<double>& A, const std::vector<double>& b, std::vector<double>& x, int N, std::vector<double>& mem)
{
    double ans = 0;

    mem = MultiplyMatrixByVector(A, x);

    for (int i = 0; i < N; ++i){
        ans += (b[i] - mem[i]) * (b[i] - mem[i]);
    }

    return sqrt(ans);
}