#include "funcs.hpp"


int main(int argc, char *argv[]) {
    try { 
        
        if (argc >= 3 && std::string(argv[1]) == "convergence"){
            std::string scheme = argv[2];
            int testCase = (argc >= 4) ? std::stoi(argv[3]) : 1; // По умолчанию case 1
            
            if (scheme == "explicit" || scheme == "implicit") {
                runConvergenceTest("convergence_results.txt", scheme, 1.0, testCase);
                std::cout << "Convergence test (" << scheme 
                          << ", case " << testCase << ") completed\n";
                return 0;
            }
        }

        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " N M scheme\n"
                      << "  N - number of time grid points\n"
                      << "  M - number of spatial grid points\n"
                      << "  scheme - 'explicit' or 'implicit'\n";
            return 1;
        }
        
        int N = std::stoi(argv[1]);
        int M = std::stoi(argv[2]);
        std::string scheme(argv[3]);

        double T = 1.; 
        double h = 1.0 / (M - 0.5);
        // double tau = h*h/4;  
        double tau = T/N;

        int levels = findSteps(T,tau);
        // if ((T - levels*tau) > 1e-6){
        //     levels++;
        // }
        N = levels;

        // std::cout<<"NUM OF STEPS = "<<levels<<std::endl;

        double* xGrid = new double[M+1];
        for (int j = 0; j < M+1; ++j) {
            xGrid[j] = -h/2. + j*h;
        }
        saveVectorToFile(xGrid, M+1, "x_grid.txt");
        delete[] xGrid;

        double* tGrid = new double[N+1];
        for (int i = 0; i < N+1; ++i) {
            tGrid[i] = i * tau;
        }
        saveVectorToFile(tGrid, N+1, "t_grid.txt");
        delete[] tGrid;

        if (N < 3 || M < 3) {
            throw std::invalid_argument("Grid sizes must be at least 3");
        }
        
        // double** u = nullptr;
        double** u = new double* [N+1];
        for (int i = 0; i < N+1; ++i) {
            u[i] = new double[M + 1];  
            for (int j = 0; j < M+1; ++j) {
                u[i][j] = -22.0;  
            }
        }

        // printMatr(u, N, M);//////////////////////

        // первая координата будет двигаться по переменной t, вторая по x
        if (scheme == "implicit") {
            solveImplicit(M, N, tau, u);
        } 
        else if (scheme == "explicit") {
            solveExplicit(N, M, tau, h, u); 
        } 
        else {
            throw std::invalid_argument("Invalid scheme, use 'explicit' or 'implicit'");
        }

        saveMatrixToFile(u, N+1, M+1, "solution_matrix.txt");

        
        std::cout<<"error = "<<error(u, M, N, h, T)<<std::endl;


        freeSolution(u, N);


    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return -1;
    }

    return 0;
}