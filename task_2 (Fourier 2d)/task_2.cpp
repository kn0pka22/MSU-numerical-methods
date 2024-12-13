#include "task_2.hpp"


double f(double x, double y){
    //return x * (1 - x) * y * (1 - y) * cos(x * x) * cos(y * y);
    //return x*(1-x)*y*(1-y);
    return(sin(3*M_PI*x)*sin(3*M_PI*y));
}


void FillingNodes(double* xk, int N){
    double h = 1/((double)N-0.5);
    xk[0] = -h/2.;
    for (int i=1;i<N+1;++i){
        xk[i] = xk[i-1]+h;
    }
}


void PrintMatrix(const double* matrix, int N, const std::string& name) {
    std::cout << "==================================================================" << name << "==================================================================" << std::endl;

    for (int k = 0; k < N * N; k++) {
        if (k % N == 0 && k != 0) std::cout << std::endl;
        std::cout << std::fixed << std::setprecision(15) << std::setw(20) << *(matrix + k) << "    "; 
    }
    std::cout << std::endl << std::endl;
}

void FillingUMatrix(int N, double* U, double* xk, double (*f)(double, double)){
    for (int i = 0; i < N+1; ++i){
        for (int j = 0; j < N+1; ++j){
            U[i * (N+1) + j] = f(xk[i], xk[j]);
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

double ScalarProduct(double* ar1,double* ar2, int N){  //подумать тут ещё над точностью
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

void FillingDMatrix(int N, double* D, double* U, double* phi){
    for (int j = 1; j < N; ++j){ 
        D[j * (N+1)] = 0;
        D[j * (N+1) + N] = 0;
        CoeffCalculate(N, U + j * (N+1), phi, D + j * (N+1));
    }
    //PrintMatrix(D, N+1, "D");
}

void FillingCMatrix(int N, double* D, double* C, double* fmemory, double* phi){
    for (int i = 1; i < N; ++i){
        C[i * (N+1)] = 0;
        C[i * (N+1) + N] = 0;
        for (int h = 0; h < N+1; ++h){
            fmemory[h] = D[h *(N+1) + i];
        }
        CoeffCalculate(N, fmemory, phi, C + i * (N+1));  
    }
    //PrintMatrix(C, N+1, "C");
}



double Calc2DFourier(double* C, int N, double x, double y){
    double res = 0;
    for (int m = 1; m < N; ++m){
        for (int n = 1; n < N; ++n){
            res += C[m * (N+1) + n] * sin(M_PI * n * x) * sin(M_PI * m * y);
        }
    }
    return res;
}





double WriteToConsole(int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi){
    //double h = 1/(N-0.5);
    double xi = 0;
    double yi = 0;
    double deltax = 0;
    double deltay = 0;
    std::cout<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
    <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
    <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
    <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;

    clock_t start=clock();
    FillingNodes( xk, N);
    FillingUMatrix(N, U, xk, f);
    FillingDMatrix(N, D, U,phi);
    FillingCMatrix(N, D, C, fmemory, phi);
    clock_t end=clock();
    double duration =(double)(end-start)/CLOCKS_PER_SEC;
    

    for (int i = 1; i < N-1; ++i){ 
        for (int j = 1; j < N; ++j){ 
            xi = xk[i];
            deltax = xk[i+1] - xi;
            deltax /= 3.;    
            yi = xk[j];
            deltay = xk[j+1] - yi;
            deltay /= 3.;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

            xi += deltax;   
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

            xi += deltax;   
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;
        }
    }
    return duration;
}


void  WriteToFile(const std::string& filename, int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi){
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        //double h = 1/(N-0.5);
        double xi = 0;
        double yi = 0;
        double deltax = 0;
        double deltay = 0;
        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
        <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;
        FillingNodes( xk, N);
        FillingUMatrix(N, U, xk, f);
        FillingDMatrix(N, D, U,phi);
        FillingCMatrix(N, D, C, fmemory, phi);
        for (int i = 1; i < N-1; ++i){ 
            for (int j = 1; j < N; ++j){ 
                xi = xk[i];
                deltax = xk[i+1] - xi;
                deltax /= 3.;    
                yi = xk[j];
                deltay = xk[j+1] - yi;
                deltay /= 3.;
                outFile << std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

                xi += deltax;   
                yi += deltay;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

                xi += deltax;   
                yi += deltay;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;
            }
        }
    }
    else {
        std::cerr << "Error opening file" << std::endl;
    }
}

double normFunction(double (*f)(double, double), double *C, int N){
    double h = 1./100.;
    double max = -1.;
    double delta = 0.;
    for (double x = 0.; x < 1.; x += h){
        for (double y = 0.; y < 1.; y += h){
            delta = fabs(f(x, y) - Calc2DFourier(C,N, x, y));
            //std::cout<<delta<<std::endl;
            if (delta > max)
                max = delta;
            //if (max>0.5) break;
        }
    }
    return max;
}
