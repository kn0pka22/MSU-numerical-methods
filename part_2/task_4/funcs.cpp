#include "funcs.hpp"


double f(double t, double x, double y){
    double p = 0.;
    //return x * (1 - x) * y * (1 - y) * cos(x * x) * cos(y * y);
    //return x*(1-x)*y*(1-y);
    // return(sin(3*M_PI*x)*sin(3*M_PI*y) + 1.*t);
    // return 0*t + 0*x + 0*y;
    // return exp(-M_PI*M_PI*(1+1)*t) * sin(M_PI*x)*sin(M_PI*y);
    return exp(-M_PI * M_PI * 0.25 * t) * ((-M_PI * M_PI * 0.25) * (x*x - x) * (y*y - y) - 2. * (y*y - y) - 2. * (x*x - x) + p* ((x*x - x) * (y*y - y)));
    // return 0*x+0*y+0*t;
}

double u0(double x, double y){
    // return sin(M_PI*x)*sin(M_PI*y);
    // return 0.*x + 0.*y;
    return (x*x - x) * (y*y - y);
    // return sin(M_PI*x)*sin(2.*M_PI*y);
}

double u0WithT(double t, double x, double y){
    // return sin(M_PI*x)*sin(M_PI*y);
    // return 0.*x + 0.*y;
    return (x*x - x) * (y*y - y) + 0*t;
    // return sin(M_PI*x)*sin(2.*M_PI*y) + 0*t;
}



double uFunc(double t, double x, double y){
    // return exp(-M_PI*M_PI*(1+1)*t) * sin(M_PI*x)*sin(M_PI*y);
    return exp(-M_PI * M_PI * 0.25 * t) * ((x*x - x) * (y*y - y));
    // return sin(M_PI*x)*sin(2.*M_PI*y)*exp(-5.*M_PI*M_PI*t);
}



double uuuFunc(double* b, int i, int j, int Nx){
    // return exp(-M_PI*M_PI*(1+1)*t) * sin(M_PI*x)*sin(M_PI*y);
    // return exp(-M_PI * M_PI * 0.25 * t) * ((x*x - x) * (y*y - y));
    // return sin(M_PI*x)*sin(2.*M_PI*y)*exp(-5.*M_PI*M_PI*t);
    return b[i*(Nx+1)+j];
}


void FillingNodes(double* xk, int N){
    double h = 1/((double)N-0.5);
    xk[0] = -h/2.;
    for (int i=1;i<N+1;++i){
        xk[i] = xk[i-1]+h;
    }
}

double p(double x, double y){
    return 0*x + 0*y + 0.;
}

void FillingNodesForTime(double* tk, int Nt){
    double tau = 0.1/Nt;
    tk[0] = 0.;
    for (int i=1;i<Nt+1;++i){
        tk[i] = tk[i-1]+tau;
    }
}


void PrintMatrix(const double* matrix, int rows, int columns, const std::string& name) {
    std::cout << "==================================================================" 
              << name 
              << "==================================================================" 
              << std::endl;
    
    for (int m = 1; m < rows-1; ++m) {
        for (int n = 1; n < columns-1; ++n) {
            std::cout << std::fixed << std::setprecision(15) 
                      << std::setw(20) << matrix[m * columns + n] << "    "; 
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

void FillingFMatrix(int Nx, int Ny, double* F, double t, double* xk, double* yk, double (*f)(double, double, double)){
    for (int i = 0; i < Ny+1; ++i){
        for (int j = 0; j < Nx+1; ++j){
            F[i * (Nx+1) + j] = f(t, xk[j], yk[i]);
        }
    }
    // PrintMatrix(F, Ny+1, Nx+1, "F");
}

void FillingUMatrix(int N, double* U, double* xk, double t, double (*u)(double, double, double)) {
    for (int i = 0; i < N+1; ++i){
        for (int j = 0; j < N+1; ++j){
            U[i * (N+1) + j] = u(t, xk[i], xk[j]);
        }
    }
}


void FillingU0Matrix(int Nx, int Ny, double* U, double* xk, double* yk, double (*u0)(double, double)) {
    for (int i = 0; i < Ny+1; ++i){
        for (int j = 0; j < Nx+1; ++j){
            U[i * (Nx+1) + j] = u0(xk[j], yk[i]);
        }
    }
}


double PhiCalculate(int n, int k, double h){
    return sin(M_PI * n * (-h/2. + k*h));
}

void PhiVectorCalculate(int N, int n, double* phi){  //phi_k^{(n)} for fixed n
    double h = 1/(N-0.5);
    for (int k=1;k<N;++k){
        phi[k] = PhiCalculate(n, k, h);   
    }
}

double ScalarProduct(double* ar1,double* ar2, int N){  
    double h = 1/((double)N-0.5);
    double res = 0.;
    for (int k=1;k<N; ++k){
        res += ar1[k] * ar2[k] * h;  
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
        cn[n] = a/b; //a*2.;
    }  
}

void FillingDMatrix(int Nx, int Ny, double* D, double* U, double* phi){
    for (int i = 0; i<Nx+1; ++i){
        D[0  * (Nx+1) + i] = 0.;
        D[Ny * (Nx+1) + i] = 0.;
    }
    for (int j = 1; j < Ny; ++j){ 
        D[j * (Nx+1)] = 0;
        D[j * (Nx+1) + Nx] = 0;
        CoeffCalculate(Nx, U + j * (Nx+1), phi, D + j * (Nx+1));
    }
    // PrintMatrix(D, Ny+1, Nx+1, "D");

}

void FillingCMatrix(int Nx, int Ny, double* D, double* C, double* fMem, double* phi){

    // матрица С имеет размерность (Nx+1)x(Ny+1) (в отличие от D, которая имеет размерность  (Ny+1)х(Nx+1) )
    for (int i = 1; i < Nx; ++i){
        C[i * (Ny+1)] = 0;
        C[i * (Ny+1) + Ny] = 0;
        for (int h = 0; h < Ny+1; ++h){  // Nx+1 -- количество столбцов в матрице D
            fMem[h] = D[h *(Nx+1) + i];
        }
        CoeffCalculate(Ny, fMem, phi, C + i * (Ny+1));  
    }
    // PrintMatrix(C, Nx+1, Ny+1, "C");

    //PrintMatrix(C, N+1, "C");
}


double Calc1DFourier(double* C, int N, double x){
    double res = 0.;
    for (int m = 1; m < N; ++m){
        res += C[m] * sin(M_PI * m * x);
    }
    return res;
}

double Calc2DFourier(double* C, int Nx, int Ny, double x, double y){
    double res = 0;
    for (int m = 1; m < Nx; ++m){
        for (int n = 1; n < Ny; ++n){
            res += C[m * (Ny+1) + n] * sin(M_PI * m * x) * sin(M_PI * n * y);
        }
    }
    return res;
}

void CalculateFourierCoefficients(int Nx, int Ny, double t, double* xk, double* yk, 
                                double* level, double* D, double* C, 
                                double* phi, double* fMem, 
                                double (*f)(double, double, double)){
    FillingNodes(xk, Nx);
    FillingNodes(yk, Ny);
    FillingFMatrix(Nx, Ny, level, t, xk, yk, f);
    FillingDMatrix(Nx, Ny, D, level, phi);
    FillingCMatrix(Nx, Ny, D, C, fMem, phi);
}



void CalculateFourierCoefficientsWithReadyMatrix(int Nx, double Ny, double* matr, double* D, double* C, 
                                            double* phi, double* fMem){
    FillingDMatrix(Nx, Ny, D, matr, phi);
    FillingCMatrix(Nx, Ny, D, C, fMem, phi);
}

double lambda(int N, double n){
    double h = 1./((double)N - 0.5);
    return (4./(h*h))*sin(M_PI*n*h/2.)*sin(M_PI*n*h/2.);
}

void CalculateCoeffsForNextLevel(int Nx, int Ny, double p, double tau, double* fCoeffs, double* resCoeffs){
    for(int i = 0; i < Nx+1; ++i){
        for(int j = 0; j < Ny+1; ++j){ 
            resCoeffs[i*(Ny+1)+j] = (resCoeffs[i*(Ny+1)+j]/tau + fCoeffs[i*(Ny+1)+j]) / (lambda(Nx,i) + lambda(Ny,j) + p + 1./tau);
            // resCoeffs[i*(Ny+1)+j] = (resCoeffs[i*(Ny+1)+j]/tau + fCoeffs[i*(Ny+1)+j]) / (lambda(Nx,i) + lambda(Ny,j) + p + 1./tau);
            
        }
    }
}


// void coeffsForf(double (*f)(double, double), double* U, double* xk, int N, double *C, double tau, int layer, double* phi, double* D, double* fmemory, double* d_mn, double* U_layer){

//     // FillingNodes(xk, N);
//     FillingUMatrix(N, U, xk, f, U_layer);
//     FillingDMatrix(N, D, U, phi);
//     FillingCMatrix(N, D, C, fmemory, phi);

//     for (int m=0;m<N+1;++m){
//         for (int n=0;n<N+1;++n){
//             d_mn[m*(N+1)+n] = U[m*(N+1)+n]/tau + C[m * (N+1) + n];
//         }
//     }
// }

// void newLevelU(double (*f)(double, double), double* U, double* xk, int N, double *C, double tau, int layer, double* phi, double* D, double* fmemory, double* d_mn, double* U_layer){
//     double h = 1./(N-0.5);
//     for (int m=0;m<N+1;++m){
//         for (int n=0;n<N+1;++n){
//             U[m*(N+1)+n] = (U[m*(N+1)+n]/tau + d_mn[m * (N+1) + n])/  \
//                 (4/(h*h)*sin(M_PI*m*(-h/2.+n*h)*sin(M_PI*m*(-h/2.+n*h))))+
//                 (4/(h*h)*sin(M_PI*m*(-h/2.+n*h)*sin(M_PI*m*(-h/2.+n*h))))+
//                 1./tau ;

            
//         }
//     }
// }

// void FindFourierCoefs(double* U, double* D, double* C, int N,
//                       double* fMem, double* xk, double (*u)(double, double),
//                       double* netmemory, double* uMem, double* phi){
    
//     FillingNodes(xk, N);
//     FillingUMatrix(N, U, xk, u);
//     FillingDMatrix(N, D, U, phi);
//     FillingCMatrix(N, D, C, fMem, phi);
// }


// void solvePDE2d(double (*u0)(double, double),
//                 double (*f)(double, double, double), 
//                 int N, int Nt,
//                 double* uijn,
//                 double* U, double* D, double* Cmatrix,
//                 double* xk, double* fMem, double* phi) {
    
//     double h = 1.0/(N - 0.5);
//     double tau = 1.0 / Nt;
    
//     // FillingNodes(xk, N);
//     // FillingU0Matrix(N, U, xk, u0);
//     // FillingDMatrix(N, D, U, phi);
//     // FillingCMatrix(N, D, Cmatrix, fMem, phi);

//     // Сохраняем начальный слой (t=0)
//     for (int i = 0; i < (N+1)*(N+1); ++i) {
//         uijn[i] = Cmatrix[i];
//     }

//     double* F = new double[(N+1)*(N+1)];
//     double* C_f = new double[(N+1)*(N+1)];

//     for (int G = 1; G <= Nt; ++G) {
//         double t = G * tau;
//         // FillingFMatrix(N, F, t, xk, xk, f);
//         // FillingDMatrix(N, D, F, phi);
//         // FillingCMatrix(N, D, C_f, fMem, phi);

//         for (int m = 0; m < N+1; ++m) {
//             for (int n = 0; n < N+1; ++n) {
//                 // if (m == 0 || m == N || n == 0 || n == N) {
//                 //     uijn[G*(N+1)*(N+1) + m*(N+1)+n] = 0.0;
//                 // } else {
//                     double sin_m = sin(M_PI * m * h / 2.0);
//                     double sin_n = sin(M_PI * n * h / 2.0);
//                     double denom = (4.0/(h*h))*(sin_m*sin_m + sin_n*sin_n) + 1.0/tau;
//                     double prev_coef = uijn[(G-1)*(N+1)*(N+1) + m*(N+1)+n];
//                     uijn[G*(N+1)*(N+1) + m*(N+1)+n] = 
//                         (prev_coef / tau + C_f[m*(N+1)+n]) / denom;
//                 // }
//             }
//         }
//     }

//     // Проверка результатов на последнем временном слое
//     double t_final = Nt * tau;
//     double* C_final = uijn + Nt * (N+1)*(N+1); // Коэффициенты для последнего слоя
    
//     std::cout << "\nValidation at final time t = " << t_final << ":\n";
//     std::cout << std::setw(20) << "x" << std::setw(20) << "y" 
//               << std::setw(25) << "Exact Solution" 
//               << std::setw(25) << "Computed Solution" 
//               << std::setw(20) << "Error" << std::endl;
    
//     double max_error = 0.0;
//     for (int i = 1; i < N; i++) {
//         for (int j = 1; j < N; j++) {
//             double x_val = xk[i];
//             double y_val = xk[j];
//             double exact = uFunc(t_final, x_val, y_val);
//             double computed = Calc2DFourier(C_final, N, x_val, y_val);
//             double error = fabs(exact - computed);
            
//             if (error > max_error) max_error = error;
            
//             std::cout << std::setprecision(6) << std::fixed
//                       << std::setw(20) << x_val
//                       << std::setw(20) << y_val
//                       << std::setw(25) << exact
//                       << std::setw(25) << computed
//                       << std::setw(20) << error << std::endl;
//         }
//     }
    
//     std::cout << "\nMaximum absolute error: " << max_error << std::endl;
//     std::cout << "N = " << N << ", Nt = " << Nt 
//               << ", h = " << h << ", tau = " << tau << std::endl;

//     delete[] F;
//     delete[] C_f;
// }


// double WriteToConsole(int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi){
//     //double h = 1/(N-0.5);
//     double xi = 0;
//     double yi = 0;
//     double deltax = 0;
//     double deltay = 0;
//     std::cout<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
//     <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
//     <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
//     <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;

//     clock_t start=clock();
//     FillingNodes( xk, N);
//     FillingUMatrix(N, U, xk, f);
//     FillingDMatrix(N, D, U,phi);
//     FillingCMatrix(N, D, C, fmemory, phi);
//     clock_t end=clock();
//     double duration =(double)(end-start)/CLOCKS_PER_SEC;
    

//     for (int i = 1; i < N-1; ++i){ 
//         for (int j = 1; j < N; ++j){ 
//             xi = xk[i];
//             deltax = xk[i+1] - xi;
//             deltax /= 3.;    
//             yi = xk[j];
//             deltay = xk[j+1] - yi;
//             deltay /= 3.;
//             std::cout << std::setprecision(15) << std::fixed \
//             << std::setw(20) << xi << " " \
//             << std::setw(20) << yi << " " \
//             << std::setw(20) << f(xi,yi) << " " \
//             << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

//             xi += deltax;   
//             yi += deltay;
//             std::cout << std::setprecision(15) << std::fixed \
//             << std::setw(20) << xi << " " \
//             << std::setw(20) << yi << " " \
//             << std::setw(20) << f(xi,yi) << " " \
//             << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

//             xi += deltax;   
//             yi += deltay;
//             std::cout << std::setprecision(15) << std::fixed \
//             << std::setw(20) << xi << " " \
//             << std::setw(20) << yi << " " \
//             << std::setw(20) << f(xi,yi) << " " \
//             << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;
//         }
//     }
//     return duration;
// }

void WriteToFile(std::string flag, const std::string& filename, int Nx, int Ny, int Nt, 
                double* tk, double* xk, double* yk, double* U, double* C, double* D, 
                double* fmemory, double* phi, 
                double (*f_func)(double, double, double), bool onlyNodes) {  
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    try {
        outFile << std::setw(20) << "t" 
                << std::setw(20) << "x" 
                << std::setw(20) << "y" 
                << std::setw(20) << "f(t,x,y)" 
                << std::setw(20) << "Fourier" 
                << std::setw(20) << "Error" 
                << std::endl;

        double max_error = 0.0;
        double max_denom = 0.0;


        if (flag == "txy") {
            for (int i = 1; i < Nt-1; ++i) { 
                double ti = tk[i];
                double deltat = onlyNodes ? 0 : (tk[i+1] - ti) / 3.0;

                for (int j = 1; j < Nx-1; ++j) { 
                    double xi = xk[j];
                    double deltax = onlyNodes ? 0 : (xk[j+1] - xi) / 3.0;

                    for (int k = 1; k < Ny-1; ++k) {
                        double yi = yk[k];
                        double deltay = onlyNodes ? 0 : (yk[k+1] - yi) / 3.0;

                        int steps = onlyNodes ? 1 : 3;
                        for (int n = 0; n < steps; ++n) {
                            double exact_value = f_func(ti, xi, yi); 
                            double fourier_value = Calc2DFourier(C, Nx, Ny, xi, yi);
                            double error = std::abs(exact_value - fourier_value);
                            double denom = std::abs(exact_value);

                            
                        if (error > max_error) {
                            max_error = error;
                        }
                        if (denom > max_denom) {
                            max_denom = denom;
                        }
                        // std::cout<<"MAX DENOM =============== "<<max_denom<<std::endl;

                            outFile << std::setprecision(15) << std::fixed
                                   << std::setw(20) << ti << " "
                                   << std::setw(20) << xi << " "
                                   << std::setw(20) << yi << " "
                                   << std::setw(20) << exact_value << " "
                                   << std::setw(20) << fourier_value << " "
                                   << std::setw(20) << max_error/max_denom
                                   << std::endl;

                            ti += deltat;
                            xi += deltax;
                            yi += deltay;
                        }
                    }
                }
            }
        } else {
            // t=0
            double t = 0;
            for (int j = 1; j < Nx-1; ++j) { 
                double xi = xk[j];
                double deltax = onlyNodes ? 0 : (xk[j+1] - xi) / 3.0;

                for (int k = 1; k < Ny-1; ++k) {
                    double yi = yk[k];
                    double deltay = onlyNodes ? 0 : (yk[k+1] - yi) / 3.0;

                    int steps = onlyNodes ? 1 : 3;
                    for (int n = 0; n < steps; ++n) {
                        double exact_value = f_func(t, xi, yi);  
                        double fourier_value = Calc2DFourier(C, Nx, Ny, xi, yi);
                        double error = std::abs(exact_value - fourier_value);\
                        double denom = std::abs(exact_value);
                        
                        if (error > max_error) {
                            max_error = error;
                        }
                        if (denom > max_denom) {
                            max_denom = denom;
                        }
                        std::cout<<"MAX DENOM =============== "<<max_denom<<std::endl;

                        outFile << std::setprecision(15) << std::fixed
                               << std::setw(20) << t << " "
                               << std::setw(20) << xi << " "
                               << std::setw(20) << yi << " "
                               << std::setw(20) << exact_value << " "
                               << std::setw(20) << fourier_value << " "
                               << std::setw(20) << max_error/max_denom
                               << std::endl;

                        xi += deltax;
                        yi += deltay;
                    }
                }
            }
        }

        outFile << "\nMaximum  error: " << std::setprecision(15) << max_error/max_denom << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
    }

    outFile.close();
}

void WriteToFileSimple(const std::string& filename, int Nx, int Ny, 
                double tFinal, double* xk, double* yk, double* C, double* D, 
                double* fmemory, double* phi, 
                double (*f_func)(double, double, double)) {  
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    try {
        outFile << std::setw(20) << "t" 
                << std::setw(20) << "x" 
                << std::setw(20) << "y" 
                << std::setw(20) << "f(t,x,y)" 
                << std::setw(20) << "Fourier" 
                << std::setw(20) << "Error" 
                << std::endl;

        double max_error = 0.0;
        double max_denom = 0.0;


        for (int i=1; i<Ny; ++i){
            for (int j=1; j<Nx; ++j){
                double exact_value = f_func(tFinal, xk[j], yk[j]); 
                double fourier_value = Calc2DFourier(C, Nx, Ny, xk[j], yk[j]);
                double error = std::abs(exact_value - fourier_value);
                double denom = std::abs(exact_value);
                    if (error > max_error) {
                        max_error = error;
                    }
                    if (denom > max_denom) {
                        max_denom = denom;
                    }
                    outFile << std::setprecision(15) << std::fixed
                            << std::setw(20) << tFinal << " "
                            << std::setw(20) << xk[j] << " "
                            << std::setw(20) << yk[i] << " "
                            << std::setw(20) << exact_value << " "
                            << std::setw(20) << fourier_value << " "
                            << std::setw(20) << max_error/max_denom
                            << std::endl;

            }
        }

        outFile << "\nMaximum  error: " << std::setprecision(15) << max_error/max_denom << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
    }

    outFile.close();
}
// double normFunction(double (*f)(double, double), double *C, int N){
//     double h = 1./100.;
//     double max = -1.;
//     double delta = 0.;
//     for (double x = 0.; x < 1.; x += h){
//         for (double y = 0.; y < 1.; y += h){
//             delta = fabs(f(x, y) - Calc2DFourier(C,N, x, y));
//             //std::cout<<delta<<std::endl;
//             if (delta > max)
//                 max = delta;
//             //if (max>0.5) break;
//         }
//     }
//     return max;
// }






void BSolver(int Nx, int Ny, int Nt, double p, double theta, 
            double* level, double* xk, double* yk, double* b, double* coeffsForU, double* coeffsForF,
            double* prod, double* y, double* phi, double* fMem, double* D){

    double err = 1024.;
    double tau = 0.1/(double)Nt;
    theta = 1;

    int numOfIterations = 1;
    // double* u_prev = new double[(Ny+1)*(Nx+1)];


    while(numOfIterations < 3 && err > 1e-7){

        AMultiplyByX(Nx, Ny, tau, xk, yk, level, prod);
        // PrintMatrix(prod, Ny+1, Nx+1, "prod");

        for(int i = 0; i < Ny+1; ++i){
            for(int j = 0; j < Nx+1; ++j){
                prod[i*(Nx+1)+j] = b[i*(Nx+1)+j] - prod[i*(Nx+1)+j];
            }
        }


        CalculateFourierCoefficientsWithReadyMatrix(Nx, Ny, prod, D, coeffsForF, phi, fMem);
        // PrintMatrix(coeffsForF, Nx+1, Ny+1, "coeffsForF");

        
        for (int i = 0; i < Ny+1; ++i){
            for(int j = 0; j < Nx+1; ++j){ 
                coeffsForU[i*(Nx+1)+j] = (coeffsForF[i*(Nx+1)+j]) / (lambda(Nx, j) + lambda(Ny, i) + p + 1./tau);
            }
        }


        for (int i = 0; i < Ny+1; ++i){
            for(int j = 0; j < Nx+1; ++j){
                y[i*(Nx+1)+j] = Calc2DFourier(coeffsForU, Nx, Ny, xk[j], yk[i]);
            }
        }


        for (int i = 0; i < Ny+1; ++i){
            for(int j = 0; j < Nx+1; ++j){
                level[i*(Nx+1)+j] = level[i*(Nx+1)+j] + theta * y[i*(Nx+1)+j];
            }
        }
        


        err = error(Nx,Ny,Nt,xk,yk,level,b,prod);

        // std::cout<<"after error"<<std::endl;

        numOfIterations++;

    }
    std::cout<<"iteration =  "<<numOfIterations<<", error = "<<err<<std::endl;


    // delete[] u_prev;
 }

double error(int Nx, int Ny, int Nt, double* xk, double* yk, double* u, double* b, double* prod){

    double err = 0.;
    double hx = 1./((double)Nx-0.5);
    double hy = 1./((double)Ny-0.5);
    double tau = 0.1/(double)Nt;

    double diff;

    AMultiplyByX(Nx, Ny, tau, xk, yk, u, prod);

    for (int i = 1; i < Ny; ++i){
            for (int j = 1; j < Nx; ++j){
                diff = fabs(prod[i*(Nx+1)+j] - b[i*(Nx+1)+j]);
                if (diff > err) err = diff;
            }
        }

    return err;
}

void AMultiplyByX(int Nx, int Ny, double tau, double* xk, double* yk, double* u, double* res){

    double hx = 1./((double)Nx - 0.5);
    double hy = 1./((double)Ny - 0.5);


    // for(int i=1;i<Ny;++i){
    //     for(int j=1; j<Nx;++j){
            
              
            
            // 4 угла:
            // if ((i==1)&&(j==1))
            //     res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 3.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 3.*u[i*(Nx+1)+j])/(hy*hy);

            // else if((i == Ny-1) && (j == Nx-1))
            //    res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j-1)] - 2.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i-1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j])/(hy*hy);
            
            // else if ((i == 1) && (j == Nx-1))
            //    res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j-1)] - 2.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 3.*u[i*(Nx+1)+j])/(hy*hy);
            // //    res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j-1)] - 2.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j])/(hy*hy);

                                  
            // else if((i == Ny-1) && (j == 1))
            //    res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 3.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i-1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j])/(hy*hy);
            // //    res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i-1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j])/(hy*hy);



            // // стенки
            // else if ((i==1))
            //     res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j] + u[i*(Nx+1)+(j-1)])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 3.*u[i*(Nx+1)+j])/(hy*hy);
                                  
            // else if ((j==1))
            //     res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 3.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j] + u[(i-1)*(Nx+1)+j])/(hy*hy);
                  
            // else if ((i == Ny-1))
            //     res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j] + u[i*(Nx+1)+(j-1)])/(hx*hx))-(u[(i-1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j])/(hy*hy);
                                  
            // else if ((j == Nx-1))
            //     res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j-1)] - 2.*u[i*(Nx+1)+j])/(hx*hx))-(u[(i+1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j] + u[(i-1)*(Nx+1)+j])/(hy*hy);
            
            // else 
                // res[i*(Nx+1)+j] = -((u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j] + u[i*(Nx+1)+(j-1)])/(hx*hx))
                //                   -(u[(i+1)*(Nx+1)+j]  - 2.*u[i*(Nx+1)+j] + u[(i-1)*(Nx+1)+j])/(hy*hy)
                //                   + u[i*(Nx+1)+j]/tau 
                //                   + p(xk[j], yk[i]) * u[i*(Nx+1)+j];
                    
            // res[i*(Nx+1)+j] += u[i*(Nx+1)+j]/tau + p(xk[j], yk[i]) * u[i*(Nx+1)+j];
            
            // if ((i==0) || (j==0) || (i==Ny) || (j==Nx)) 
            //     res[i*(Nx+1)+j] = 0.;


        for(int i = 0; i < Ny+1; ++i) {
            for(int j = 0; j < Nx+1; ++j) {

                // if ((i == 0) || (j == 0) || (i == Ny) || (j == Nx)) {
                //     res[i*(Nx+1)+j] = 0.0;
                //     continue;
                // }
                double laplacian_x = 0.0;
                double laplacian_y = 0.0;   
                            
                if (j == 1) {
                    // Левая граница: u(i,0) = -u(i,1)
                    laplacian_x = (u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j] - u[i*(Nx+1)+j])/(hx*hx);
                } 
                else if (j == Nx-1) {
                    // Правая граница: u(i,Nx) = 0
                    laplacian_x = (0. - 2.*u[i*(Nx+1)+j] + u[i*(Nx+1)+(j-1)])/(hx*hx);
                } 
                else {
                    // Внутренние точки по x
                    laplacian_x = (u[i*(Nx+1)+(j+1)] - 2.*u[i*(Nx+1)+j] + u[i*(Nx+1)+(j-1)])/(hx*hx);
                }

                // По y
                if (i == 1) {
                    // Нижняя граница: u(0,j) = -u(1,j)
                    laplacian_y = (u[(i+1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j] - u[i*(Nx+1)+j])/(hy*hy);
                } 
                else if (i == Ny-1) {
                    // Верхняя граница: u(Ny,j) = 0
                    laplacian_y = (0. - 2.*u[i*(Nx+1)+j] + u[(i-1)*(Nx+1)+j])/(hy*hy);
                } 
                else {
                    // Внутренние точки по y
                    laplacian_y = (u[(i+1)*(Nx+1)+j] - 2.*u[i*(Nx+1)+j] + u[(i-1)*(Nx+1)+j])/(hy*hy);
                }

                res[i*(Nx+1)+j] = - laplacian_x - laplacian_y + u[i*(Nx+1)+j] * (1./tau + p(xk[j], yk[i]));
            }
        }
    //     }
    // }
}


void WriteSolutionWithError(int Nx, int Ny, double t, double* xk, double* yk, 
                          double* level, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(6);
    file << "Слой " << t << std::endl;
    
    file << "Численное решение:" << std::endl;
    for (int i = 1; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            file << std::setw(12) << level[i*(Nx+1)+j];
        }
        file << std::endl; 
    }
    file << std::endl;
    
    file << "Точное решение:" << std::endl;
    for (int i = 1; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            file << std::setw(12) << uFunc(t, xk[j], yk[i]);
        }
        file << std::endl; 
    }
    file << std::endl;
    
    file << "Относительная ошибка:" << std::endl;
    double max_error = 0.0;
    double max_value = 0.0;
    
    for (int i = 1; i < Ny-1; ++i) {
        for (int j = 1; j < Nx-1; ++j) {
            double exact = uFunc(t, xk[j], yk[i]);
            double numerical = level[i*(Nx+1)+j];
            double error = std::abs(exact - numerical);
            double relative_error = (exact != 0) ? error / std::abs(exact) : 0.0;
            
            file << std::setw(12) << relative_error;
            
            if (relative_error > max_error) max_error = relative_error;
            if (std::abs(exact) > max_value) max_value = std::abs(exact);
        }
        file << std::endl; 
    }
    
    file << "Максимальная относительная ошибка: " << max_error << std::endl;
    
    file << std::endl << std::endl;
    file.close();
}