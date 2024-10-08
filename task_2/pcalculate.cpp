#include "task_2.hpp"

int pcalculate(int N){

    int NumTests = 5;
    int M = (NumTests+1)*N;
    double* xk      = new double[M+1];
    double* fmemory = new double[M+1];
    double* U       = new double[(M+1)*(M+1)];
    double* D       = new double[(M+1)*(M+1)];
    double* C       = new double[(M+1)*(M+1)];
    double* phi     = new double[M+1];

    double* lognorms =  new double[NumTests+1];
    double* loghs    =  new double[NumTests+1];


 

    if ((!xk) || (!U) || (!D) || (!C) || (!phi) || (!fmemory) || (!lognorms) || (!loghs)){
        std::cout<<"Not enough memory!\n";
	 	return -1;
    }

    std::ofstream outFile("PDependence.txt");
    for (int i=0; i<NumTests; i+=2){
        double hForNorm = 1 / ((double)(2*N));
        double MaxForNorm = 1e-10;
        double ans = 0;
        double temp = 0;
        for (double x = 0; x < 1; x += hForNorm){
            for (double y = 0; y < 1; y += hForNorm){
                FillingNodes( xk, N*(i+1));
                FillingUMatrix(N*(i+1), U, xk, f);
                FillingDMatrix(N*(i+1), D, U,phi);
                FillingCMatrix(N*(i+1), D, C, fmemory, phi);
                ans = Calc2DFourier(C, N*(i+1), x, y);
                temp = fabs(ans - f(x, y));
                if (temp > MaxForNorm)
                    MaxForNorm = temp;
            }
        }
        lognorms[i] = log(1 / MaxForNorm);
        loghs[i] = log(N*(i+1));

        MaxForNorm = 1e-10;
        for (double x = 0; x < 1; x += hForNorm){
            for (double y = 0; y < 1; y += hForNorm){
                FillingNodes( xk, N*(i+2));
                FillingUMatrix(N*(i+2), U, xk, f);
                FillingDMatrix(N*(i+2), D, U,phi);
                FillingCMatrix(N*(i+2), D, C, fmemory, phi);
                ans = Calc2DFourier(C, N*(i+2), x, y);
                temp = fabs(ans - f(x, y));
                if (temp > MaxForNorm)
                    MaxForNorm = temp;
            }
        }
        lognorms[i+1] = log(1 / MaxForNorm);
        loghs[i+1] = log(N*(i+2));

        
        if (outFile.is_open()){
            outFile<<std::setw(10)<<loghs[i]<<std::setw(10)<<lognorms[i]<<std::endl;
        }
    }
    // outFile.close();

    

    double p;
    if (fabs(loghs[0] - loghs[1]) > 1e-10){
        p = (lognorms[1] - lognorms[0]) / (loghs[1] - loghs[0]);
    }
    std::cout<<" P = "<<p<<std::endl;


    delete[] xk;
    delete[] fmemory;
    delete[] U;
    delete[] D;
    delete[] C;
    delete[] phi;
    delete[] lognorms;
    delete[] loghs;

	return 0;
}