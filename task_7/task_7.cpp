#include "task_7.hpp"

double f(double x) {
    return sin( M_PI * x);
    //return x*(x-1);
    
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

void FillingNodes(std::vector<double>& xk){
    int N = xk.size(); //N+1
    N--;
    double h = 1/(N-0.5);

    xk[0] = -h/2.;
    for (int i=1;i<N+1;++i){
        xk[i] = xk[i-1]+h;
    }
}

void FillingBasicNodes(std::vector<double>& xk){
    int N = xk.size(); //N-1
    N++;
    double h = 1/(N-0.5);

    xk[0] = h/2.;
    // X_1, .., X_N-1 
    for (int i=1;i<N-1;++i){
        xk[i] = xk[i-1]+h;
    }
}


void FillingValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f){
    int N = xk.size(); //N+1
    N--;

    if (xk.empty()) {
        std::cerr<<"Please, fill nodes xk!"<<std::endl; 
        throw std::runtime_error("xk not found");
    }

    for (int i=0;i<N+1;++i){
        yk[i] = f(xk[i]);
    }
}

void FillingBasicValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f){
    int N = xk.size(); //N-1
    N++;

    if (xk.empty()) {
        std::cerr<<"Please, fill nodes xk!"<<std::endl; 
        throw std::runtime_error("xk not found");
    }

    for (int i=0;i<N-1;++i){
        yk[i] = f(xk[i]);
    }
}

double Psi(int k, int n, int N){
    return sin(M_PI * n * (k-1./2.) / (double)(N-0.5));
}

double Lambda(int n, int N, double p){
    double lam = p - 2. * (double)(N-0.5) * (double)(N-0.5) * (cos(M_PI * n / (double)(N-0.5)) - 1.);
    return lam;
}

double ScalarProduct(int n, std::vector<double>& fk){
    int N = fk.size(); //N-1
    N++;
    double sp = 0.;
    for (int i = 0; i < N-1; ++i){
        sp += 2 * fk[i] * Psi(i+1, n, N) / (double)(N-0.5);
    }
    return sp;
}


void Fourier(std::vector<double>& y, double p, std::vector<double>& fk){
    int N = fk.size(); //N-1
    N++;

    for (int k = 0; k < N - 1; ++k){
        y[k]=0.;
        for (int n = 1; n < N; ++n){
            //std::cout<<(ScalarProduct(n, fk) / Lambda(n, N, p) )* Psi(k+1, n, N)<<" ";
            y[k] += (ScalarProduct(n, fk) / Lambda(n, N, p) )* Psi(k+1, n, N);
        }
        //std::cout<<"y["<<k<<"] = "<<y[k]<<std::endl;
    }
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

void BasisMatrixFill(double p, std::vector<double>& M){
    int N = sqrt(M.size());  //N-1
    N++;

    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N-1; ++j) {
            if (i == j){
                M[i * (N - 1) + j] = p + 2.0 * (double)((N-0.5) * (N-0.5));
            } 
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N - 1) + j] = -(double)((N-0.5) * (N-0.5));
            } 
            else{
                M[i * (N - 1) + j] = 0.0;
            }
        }
    }
    M[0] = 3.* (double)((N-0.5) * (N-0.5)) + p;
}

void FullMatrixFill(double p, std::vector<double>& M){
    int N = sqrt(M.size());  //N+1
    N--;

    for (int i = 0; i < N+1; ++i) {
        for (int j = 0; j < N+1; ++j) {
            if (i == j){
                M[i * (N + 1) + j] = p + 2.0 * (double)((N-0.5) * (N-0.5));
            } 
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N + 1) + j] = -(double)((N-0.5) * (N-0.5));
            } 
            else{
                M[i * (N + 1) + j] = 0.0;
            }
        }
    }
    M[1*(N+1)+1]= p + 3.0 * (double)((N-0.5) * (N-0.5));
    M[0*(N+1)+0]= 1./2.;
    M[0*(N+1)+1]= 1./2.;
    M[N*(N+1)+N]= 1.;
    M[1*(N+1)+0]= 0.;
    M[N*(N+1)+N-1]= 0.;
    M[(N-1)*(N+1)+N]= 0.;

}


double Richardson(std::vector<double>& x, const std::vector<double>& A, const std::vector<double>& b, double tau, int n, int mIter, std::vector<double>& mem){
    for (int k = 0; k < n; ++k){
        x[k] = 0;
        mem[k] = 0;
    }
    for (int m = 0; m < mIter; m++){
        mem = MultiplyMatrixByVector(A, x);
        for (int i = 0; i < n; i++){
            x[i] = x[i] - tau * mem[i] + tau * b[i]; //-=
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


void WriteResultsToFile(const std::string& filename, std::vector<double>& x, 
                        std::vector<double>& A, const std::vector<double>& b,
                        double tau, int N, int m, double dq, double q, std::vector<double>& mem) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    file << std::left 
         << std::setw(15) << "n" 
         << std::setw(20) << "Eps" 
         << std::setw(20) << "dq*fnorm" 
         << std::endl;

    double fnorm = Richardson(x, A, b, tau, N-1, 0, mem);
    double eps = 0.0;

    for (int n = 0; n < m; ++n) {
        eps = Richardson(x, A, b, tau, N-1, n, mem);
        file << std::left 
             << std::setw(15) << std::fixed << std::setprecision(10) << n
             << std::setw(20) << std::fixed << std::setprecision(10) << eps
             << std::setw(20) << std::fixed << std::setprecision(10) << dq * fnorm
             << std::endl;
        dq *= q;
    }
}


double BSolver( std::vector<double>& x, std::vector<double>& A, 
                std::vector<double>& B, std::vector<double>& b, 
                double tau, int mIter, std::vector<double>& mem, 
                std::vector<double>& mem1, double p) {

    int N = x.size(); //N-1
    N++;

    for (int k = 0; k < N-2; ++k){
        x[k] = 0;
        mem[k] = 0;
    }
    for (int m = 0; m < mIter; ++m) {
        mem = MultiplyMatrixByVector(A, x);
        
        // b - Ax
        for (int j = 0; j < N-2; ++j) {
            mem[j] = b[j] - mem[j];
        }
        Fourier(mem1, p, mem);
        //mem1 = MultiplyMatrixByVector(B, mem);

        for (int i = 0; i < N-2; ++i) {
            x[i] += tau * mem1[i];
        }
    }

    return ErNormInf(A, b, x, mem);
}


double ErNormInf(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, std::vector<double>& mem){

    int N = x.size(); //N-1
    N++;
    double ans = 0;
    

    mem = MultiplyMatrixByVector(A, x);

    for (int i = 0; i < N-2; i++){
        if(fabs((b[i] - mem[i])) > ans) ans = fabs((b[i] - mem[i]));
    }

    return sqrt(ans);
}


double BasisMatrixFillWithVariableP(std::vector<double>& M) {
    int N = sqrt(M.size()); //N-1
    N++;
    double mean = 0;
    double pk;
    for (int i = 0; i < N-1; ++i) {
        //pk = sin(M_PI * (i - 1./2.) / (N - 0.5));
        pk = sin(M_PI * i / (N - 0.5));
        for (int j = 0; j < N-1; ++j) {
            if (i == j){
                M[i * (N - 1) + j] = 2.0 * (double)((N-0.5) * (N-0.5)) + 1.+pk*pk;
            } 
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N - 1) + j] = -(double)((N-0.5) * (N-0.5));
            } 
            else{
                M[i * (N - 1) + j] = 0.0;
            }
        }
        mean += pk*pk+1.;
    }
    mean/=(N-2);
    pk = sin(M_PI * 0 / (N - 0.5));
    M[0] = 3.* (double)((N-0.5) * (N-0.5)) + 1.+pk*pk;
    return mean;

}

double SearchQ(std::vector<double>& A){
    int N = sqrt(A.size()); //N-1
    N++;

    // double max = 0.;
    // double min = 0.;
    // double sum = 0.;
    // for (int i = 0; i < N-2; ++i){
    //     sum = 0.;
    //     for (int j = 0; j < N-2; ++j){
    //         sum += fabs(A[i * (N-2) + j]);
    //     }
    //     if(i == 0) min = 2 * A[i * (N-2) + i] - sum;
        
    //     if (sum > max) max = sum;

    //     if (2 * A[i * (N-2) + i] - sum < min) min = 2 * A[i * (N-2) + i] - sum;
    // }
    // return (double)(max - min) / (double)(max + min);

    double sum = 0;
    double qMax = 0;
    double q;
    for (int i = 0; i < N-1; ++i){
        sum = 0.;
        for (int j = 0; j < N-1; j++){
            if(j != i) sum += fabs(A[i * (N-1) + j]);;
        }
        //printMatrix(A);
        //std::cout<<sum<<" "<<fabs(A[i * (N-1) + i])<<std::endl<<std::endl;
        q = sum / fabs(A[i * (N-1) + i]);
        // if(q > 1){
        //     std::cout<<"bad q, would haven't conv\n"<<std::endl;
        //     return 0;
        // }
        if(i == 0){
            qMax = q;
        }
        if (q > qMax){
            qMax = q;
        }
    }
    return qMax;
}
