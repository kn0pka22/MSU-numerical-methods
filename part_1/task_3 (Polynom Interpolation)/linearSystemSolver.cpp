#include "linearSystemSolver.hpp"

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

int solve(std::vector<double>& M, std::vector<double>& b, std::vector<double>& x, std::vector<int>& memory) {
    int n = std::sqrt(M.size()); 

    if (n * n != M.size()) {
        std::cerr << "Error: Matrix size is not valid." << std::endl;
        return -1;
    }
    memory.resize(n);
    for (int j = 0; j < n; j++) {
        memory[j] = j;
    }

    double mean = 1e-5;
    for (int i = 0; i < n * n; ++i) {
        mean += std::fabs(M[i]); 
    }

    double eps = mean * 1e-10 / (n * n);

    for (int step = 1; step <= n; step++) {
        int indMax = step;
        double max = M[e(step, step, n)];

    
        for (int i = indMax + 1; i <= n; ++i) {
            if (std::fabs(M[e(step, i, n)]) > std::fabs(max)){
                indMax = i;
                max = M[e(step, i, n)];
            }
        }

        if (std::fabs(max)<std::fabs(eps)) {
            return -1; 
        }

        if (indMax != step) {
            for (int i = 1; i <= n; ++i) {
                std::swap(M[e(i, step, n)], M[e(i, indMax, n)]);
            }
            std::swap(memory[indMax - 1], memory[step - 1]);
        }

        for (int i = step + 1; i <= n; ++i) {
            double v = (M[e(i, step, n)] / max);
            for (int j = step; j <= n; ++j) {
                M[e(i, j, n)] -= v * M[e(step, j, n)];
            }
            b[i - 1] -= v * b[step - 1];
        }
    }

    for (int step = 1; step <= n; ++step) {
        for (int i = 1; i <= n - step; ++i) {
            b[i - 1] -= b[n - step] * M[e(i, n - step + 1, n)] / M[e(n - step + 1, n - step + 1, n)];
            M[e(i, n - step + 1, n)] = 0; 
        }
        b[n - step] /= M[e(n - step + 1, n - step + 1, n)];
        M[e(n - step + 1, n - step + 1, n)] = 1; 
    }

    for (int u = 0; u < n; ++u) {
        x[memory[u]] = b[u];
    }

    return 0; 
}

