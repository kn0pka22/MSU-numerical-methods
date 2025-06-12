#include "funcs.hpp"
// #include "BSolver.hpp"


int main(int argc, char *argv[]){

    int Nx = 5;
    int Ny = 5;
    int Nt = 3;

    int N = std::max(Nx, Ny);

    double* level = new double[(Ny+1)*(Nx+1)];
    double* F = new double[(Ny+1) * (Nx+1)];
    double* FF = new double[(Ny+1) * (Nx+1)];
    double* FFExact = new double[(Ny+1) * (Nx+1)];




    double* xk        = new double[Nx+1];
    double* yk        = new double[Ny+1];
    double* tk        = new double[Nt+1];

    double* D         = new double[(Ny+1)*(Nx+1)];

    double* C         = new double[(Nx+1)*(Ny+1)];

    double* phi       = new double[N+1];
    double* fMem      = new double[Ny+1];

    double* coeffsForU  = new double[(Nx+1)*(Ny+1)];
    double* coeffsForF  = new double[(Nx+1)*(Ny+1)];
    double* coeffsForNextLevel  = new double[(Ny+1)*(Nx+1)];

    double* b = new double[(Nx+1) * (Ny+1)];
    double* prod = new double[(Nx+1) * (Ny+1)];
    double* y = new double[(Nx+1) * (Ny+1)];



   

    double p = 0.;
    double theta = 1;

    double hx  = 1./((double)Nx-0.5);
    double hy  = 1./((double)Ny-0.5);
    double tau = 0.1/((double)Nt);

    

    

    FillingNodesForTime(tk, Nt);

    // for (int i=0;i<Nt+1;++i){
    //     std::cout<<tk[i]<<"  ";
    // }
    // std::cout<<std::endl;


    // for(int i=0;i<Ny+1;++i){    
    //     for(int j=0;j<Nx+1;++j){ 
    //         level[i*(Nx+1)+j] = u0(xk[j], yk[i]);
    //     }
    // }

    // FillingU0Matrix(Nx, Ny, level, xk, yk, u0);
    // PrintMatrix(level, Ny+1, Nx+1, "level");

    // FillingDMatrix(Nx, Ny, D, level, phi);
    //-------------------------------------------------------------- <- FillingDMatrix -- works correct!
    // std::cout<<"CHECK FOR CORRECTION WORK OF FillingDMatrix\n";
    // for (int i = 1; i<Ny; ++i){
    //     for (int k = 1; k < Nx; ++k){
    //         std::cout<<Calc1DFourier(D+i*(Nx+1), Nx, xk[k])<<"  "<<level[i*(Nx+1)+k]<<std::endl;
    //     }
    // }
    //----------------------------------------------------------------
    // FillingCMatrix(Nx, Ny, D, C, fMem, phi);
    //-------------------------------------------------------------- <- FillingСMatrix -- works correct!
    // std::cout << "CHECK FOR CORRECTION WORK OF FillingCMatrix\n";
    // for (int i = 1; i < Nx; ++i) {
    //     for (int k = 1; k < Ny; ++k) {
    //         std::cout << Calc1DFourier(C + i*(Ny+1), Ny, yk[k]) 
    //                 << "  " 
    //                 << D[k*(Nx+1) + i] << std::endl;
    //     }
    // }
    //----------------------------------------------------------------

    CalculateFourierCoefficients(Nx, Ny, tk[0], xk, yk, level, D, coeffsForU, phi, fMem, u0WithT);

    // for (int i=0; i<Ny+1; ++i){
    //         for (int j=0; j<Nx+1;++j){
    //             FF[i*(Nx+1)+j] = Calc2DFourier(coeffsForU, Nx, Ny, xk[j], yk[i]);
    //         }
    //     }

    // for (int i=0; i<Ny+1; ++i){
    //         for (int j=0; j<Nx+1;++j){
    //             FFExact[i*(Nx+1)+j] = u0WithT(0, xk[j], yk[i]);
    //         }
    //     }

    // PrintMatrix(FFExact, Ny+1, Nx+1, "f");
    // PrintMatrix(FF, Ny+1, Nx+1, "Fourier");




    // WriteToFile("txy", "out.txt", Nx, Ny, Nt, tk, xk, yk, level, coeffsForU, D, fMem, phi, u0WithT);  //<- correct
    // WriteToFileSimple("out.txt", Nx, Ny, Nt, tk, xk, yk, level, coeffsForU, D, fMem, phi, u0WithT);
    // double* tmpMatr = new double[(Ny+1)*(Nx+1)];
    // for (int i=0; i<Ny+1; ++i){
    //         for (int j=0; j<Nx+1;++j){
    //             tmpMatr[i*(Nx+1)+j] = level[i*(Nx+1)+j]/tau + f(0, xk[j], yk[i]);
    //         }
    //     }
    // // PrintMatrix(tmpMatr, Ny+1, Nx+1, "tmpMatr");
    // delete[] tmpMatr;


    // for (int i=1; i<Ny; ++i){
    //         for (int j=1; j<Nx; ++j){
    //             FFExact[i*(Nx+1)+j] = uFunc(tau, xk[j], yk[i]);
    //         }
    //     }

    // AMultiplyByX(Nx, Ny, tau, xk, yk, FFExact, prod);
    // PrintMatrix(prod, Ny+1, Nx+1, "prod");

    // for (int i=1; i<Ny; ++i){
    //     for (int j=1; j<Nx;++j){
    //         FFExact[i*(Nx+1)+j] = u0WithT(0, xk[j], yk[i])/tau + f(tau, xk[j], yk[i]);
    //     }
    // }
    // PrintMatrix(FFExact, Ny+1, Nx+1, "u0WithT(0, xk[j], yk[i])/tau + f(tau, xk[j], yk[i])");

      // for (int i=1; i<Ny; ++i){
    //         for (int j=1; j<Nx; ++j){
    //             FFExact[i*(Nx+1)+j] = uFunc(tau, xk[j], yk[i]);
    //         }
    //     }

    // AMultiplyByX(Nx, Ny, tau, xk, yk, FFExact, prod);
    // PrintMatrix(prod, Ny+1, Nx+1, "prod");

    

    
    for (int t=1; t<Nt; ++t){
        // std::cout<<"process started"<<std::endl;
        CalculateFourierCoefficients(Nx, Ny, tk[t], xk, yk, F, D, coeffsForF, phi, fMem, f);
            
        // FillingFMatrix(Nx, Ny, FF, 1*tau, xk, yk, f);

        // PrintMatrix(FF, Ny+1, Nx+1, "f");
        // PrintMatrix(F, Ny+1, Nx+1, "Fourier");

        // WriteToFileSimple("out.txt", Nx, Ny, Nt, tk, xk, yk, F, coeffsForF, D, fMem, phi, f);

        // WriteToFile("txy", "out.txt", Nx, Ny, Nt, tk, xk, yk, F, coeffsForF, D, fMem, phi, f); // <- correct for tk[t]

        // for (int i=0; i<Ny+1; ++i){
        //     for (int j=0; j<Nx+1;++j){
        //         FF[i*(Nx+1)+j] = Calc2DFourier(coeffsForF, Nx, Ny, xk[j], yk[i]);
        //     }
        // }

        // for (int i=0; i<Ny+1; ++i){
        //         for (int j=0; j<Nx+1;++j){
        //             FFExact[i*(Nx+1)+j] = f(tau, xk[j], yk[i]);
        //         }
        //     }

        // PrintMatrix(FFExact, Ny+1, Nx+1, "f");
        // PrintMatrix(FF, Ny+1, Nx+1, "Fourier");




    


        for(int i = 0; i < Ny+1; ++i){
            for(int j = 0; j < Nx+1; ++j){ 
                b[i*(Nx+1)+j] = level[i*(Nx+1)+j]/tau + f(t*tau, xk[j], yk[i]); //  ТЫ УМНИЦА!!!! ты умница!!!
            }
            for(int j = 0; j < Nx+1; ++j){ 
                level[i*(Nx+1)+j] = Calc2DFourier(coeffsForU, Nx, Ny, xk[j], yk[i]);
            }
        } 


        // -----------------------------

        // CalculateFourierCoefficientsWithReadyMatrix(Nx, Ny, b, D, C, phi, fMem);


        
        // for (int i=0; i<Ny+1; ++i){
        //     for (int j=0; j<Nx+1;++j){
        //         FF[i*(Nx+1)+j] = Calc2DFourier(C, Nx, Ny, xk[j], yk[i]);
        //     }
        // }

        // for (int i=0; i<Ny+1; ++i){
        //         for (int j=0; j<Nx+1;++j){
        //             FFExact[i*(Nx+1)+j] = b[i*(Nx+1)+j];
        //         }
        // }

        // PrintMatrix(FFExact, Ny+1, Nx+1, "b");
        // PrintMatrix(FF, Ny+1, Nx+1, "Fourier");

        // -----------------------------
        CalculateCoeffsForNextLevel(Nx, Ny, p, tau, coeffsForF, coeffsForU);

        for (int i=0; i<Ny+1; ++i){
            for (int j=0; j<Nx+1;++j){
                FF[i*(Nx+1)+j] = Calc2DFourier(coeffsForU, Nx, Ny, xk[j], yk[i]);
            }
        }

        AMultiplyByX(Nx, Ny, tau, xk, yk, FF, prod);
            for (int i=0; i<Ny+1; ++i){
                for (int j=0; j<Nx+1;++j){
                    FFExact[i*(Nx+1)+j] = b[i*(Nx+1)+j];
                }
        }
    // PrintMatrix(prod, Ny+1, Nx+1, "prod");

        PrintMatrix(prod, Ny+1, Nx+1, "prod");
        PrintMatrix(b, Ny+1, Nx+1, "b");





        

        // AMultiplyByX(Nx, Ny, tau, xk, yk, FFExact, prod);
        // PrintMatrix(prod, Ny+1, Nx+1, "prod");

        // for (int i=0; i<Ny+1; ++i){
        //     for (int j=0; j<Nx+1;++j){
        //         FFExact[i*(Nx+1)+j] = u0WithT(0, xk[j], yk[i])/tau + f(tau, xk[j], yk[i]);
        //     }
        // }

        // PrintMatrix(FFExact, Ny+1, Nx+1, "u0WithT(0, xk[j], yk[i])/tau + f(tau, xk[j], yk[i])");


        // WriteToFile("txy", "out.txt", Nx, Ny, Nt, tk, xk, yk, F, coeffsForU, D, fMem, phi, uFunc); // 

        // FillingFMatrix(Nx, Ny, F, t*tau, xk, yk, uFunc);
        // PrintMatrix(F, Ny+1, Nx+1, "F");
        // PrintMatrix(level, Ny+1, Nx+1, "level");


        
        // BSolver(Nx, Ny, Nt, p, theta, level, xk, yk, b, coeffsForU, coeffsForF, prod, y, phi, fMem, D);
        // CalculateFourierCoefficients(Nx, Ny, tk[0], xk, yk, level, D, coeffsForU, phi, fMem, u0WithT);

        // std::ofstream file("res.txt");
        // file << std::fixed << std::setprecision(6);
        // file << "Слой " << t << std::endl;
        // file << "Численное решение:" << std::endl;
        // for (int i = 1; i < Ny-1; ++i) {
        //     for (int j = 1; j < Nx-1; ++j) {
        //         file << std::setw(12) << level[i*(Nx+1)+j];
        //     }
        //     file << std::endl; 
        // }
        // file << std::endl;
        // file << "Точное решение:" << std::endl;
        // for (int i = 1; i < Ny-1; ++i) {
        //     for (int j = 1; j < Nx-1; ++j) {
        //         file << std::setw(12) << uFunc(tk[t], xk[j], yk[i]);
        //     }
        //     file << std::endl; 
        // }
        // file << std::endl << std::endl;
        // break;

    }

    // без BS
    WriteToFileSimple("out.txt", Nx, Ny, tk[Nt], xk, yk, coeffsForU, D, fMem, phi, uFunc);
    
    //  BS
    WriteSolutionWithError(Nx, Ny, tk[Nt], xk, yk, level, "outBS.txt");


    delete[] xk;
    delete[] D;
    delete[] C;
    delete[] phi;
    delete[] fMem;

    delete[] level;
    delete[] F;
    delete[] FF;
    delete[] FFExact;

    delete[] b;
    delete[] prod;
    delete[] y;

    delete[] coeffsForU;
    delete[] coeffsForF;
    delete[] coeffsForNextLevel;

    return 0;
}