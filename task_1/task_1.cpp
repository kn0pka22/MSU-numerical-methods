#include "task_1.hpp"


double f(double x){
    return((exp(x)-exp(1.0))*sin(x));
}


void FillingNodes(double* xk, int N){
    double h = 1/(N-0.5);

    xk[0] = -h/2.;
    for (int i=1;i<N;++i){
        xk[i] = xk[i-1]+h;
    }
}

void FillingValues(double* xk, double* yk, double (*f)(double), int N){

    if (!xk) { 
        std::cerr<<"Please, fill nodes xk!"<<std::endl; 
        throw std::runtime_error("xk not found");
    }

    for (int i=0;i<N+1;++i){
        yk[i] = f(xk[i]);
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

void CoeffCalculate(int N, int k, double* yk, double* phi, double* cn){
    double a = 0;
    double b = 0; 
    for (int n=1; n<N; ++n){
        PhiVectorCalculate(N, n, phi);
        a = ScalarProduct(yk, phi, N);
        b = ScalarProduct(phi,phi, N);
        cn[n] = a/b;
        
    }
}

double scalar(double *x, double*y, unsigned int N){
    double res = 0.;
    for (unsigned int i=1;i<N;++i){
        res += x[i]*y[i];
    }
    return res;
}

double FourierCompute(double* cn, int N, double x){
    double res = 0;
    for(int n = 1; n < N; ++n){
        res += cn[n] * sin( M_PI * n * x );
    }
    return res;
}


void WriteToConsole(int N, double* xk, double* yk, double* cn, double* phi){
    double h = 1/(N-0.5);
    std::cout<< "         xk                   yk                      yk*              "<<std::endl;
    for (int i = 1; i < N+1; ++i){
        CoeffCalculate(N, i, yk, phi, cn);
        FourierCompute(cn, N, (- h/2. + i* h));
        std::cout << std::setprecision(15) << std::fixed \
        << std::setw(20) << xk[i] << " " \
        << std::setw(20) << yk[i] << " " \
        << std::setw(20) << FourierCompute(cn, N, xk[i]) << std::endl;
        
    }
}

void WriteToFile(const std::string& filename, int N, double* xk, double* yk, double* cn, double* phi) {
    double h = 1/(N-0.5);
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        outFile << "         xk                   yk                      yk*              "<<std::endl;
        for (int i = 1; i < N+1; ++i){
            CoeffCalculate(N, i, yk, phi, cn);
            FourierCompute(cn, N, (- h/2. + i* h));
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xk[i] << " " \
            << std::setw(20) << yk[i] << " " \
            << std::setw(20) << FourierCompute(cn, N, xk[i]) << std::endl;
        }
        outFile.close();
        std::cout << "Information successfully written to file!" << std::endl;
    } 
    else {
        std::cerr << "Error opening file" << std::endl;
    }
}


//////////////////////////////////////////////////////////////////////

double NormFunction(double (*f)(double), double* cn, int N){
    double h = 1/10000.;
    double max = -1;
    double delt = 0;
    for (double x = 0.; x < 1.; x += h){
        delt = fabs(f(x) - FourierCompute(cn,N,x));
        //std::cout<<delt<<std::endl;
        if (delt > max)
            max = delt;
    }
    return max;
}


