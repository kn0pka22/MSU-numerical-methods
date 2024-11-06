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
    //int n=10;
    //std::cout<<"n = "<<n<<std::endl;

    if (n * n != M.size()) {
        std::cerr << "Error: Matrix size is not valid." << std::endl;
        return -1;
    }
    memory.resize(n);
    for (int j = 0; j < n; j++) {
        memory[j] = j;
    }

    double mean = 1e-2;
    for (int i = 0; i < n * n; ++i) {
        mean += std::fabs(M[i]); 
    }
    //std::cout<<"MEAN = "<<mean<<std::endl;

    double eps = mean * 1e-10 / (n * n);

    for (int step = 1; step <= n; step++) {
        int indMax = step;
        double max = M[e(step, step, n)];
        //std::cout<<"max = "<<max<<std::endl;

    
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



// int Michael_Jordan(int n, std::vector<double>& matrica, std::vector<double>& part_b, std::vector<double>& result)
// {
//     int i, leader;
//     double t,z;
//     z=norma_matrica(matrica,n);

//     for(i = 0; i < n; i++)
//     {
//         leader = lead(matrica, n, i, z);
//         if(leader == -1)
//         {
//             return -1;
//         }

//         if(i != leader)
//         {
//             for (int j=i; j<n;j++)
//             {
//                 t=matrica[i*n+j];
//                 matrica[i*n+j]=matrica[leader*n+j];
//                 matrica[leader*n+j]=t;
//             }
//             t=part_b[i];
//             part_b[i]=part_b[leader];
//             part_b[leader]=t;
//         }
        
//         delenie(matrica, part_b, n, i);
//         vychitanie(matrica, part_b, n, i);
        
//     }

//     for (int i=0;i<n;i++)
//     {
//         result[i]=part_b[i];
//     }

//     return 0;
// }

// int lead(std::vector<double>& matrica, int n, int j, double z)
// {
//     int i, leader;
//     double max;
//     leader=-1;
//     max=0.;

//     for(i = j; i < n; i++)
//     {
//         {  
//             if((fabs(A(i, j)) > max)&&(fabs(A(i,j))>EPS*z))
//             {
//                 //printf("A(j,i)=%.10e     MAX=%.10e\n",fabs(A(i,j)),z*EPS);
//                 leader = i;
//                 max = fabs(A(i, j));
//             }
//         }
//     }

//     return leader;
// }

// void delenie(std::vector<double>& matrica, std::vector<double>& part_b, int n, int i)
// {
//     int j;
//     double kef;
//     kef=A(i, i);
//     for(j = 0; j < n; j++)
//     {
//         A(i, j) = A(i, j)/ kef;
//     }
//     B(i)=B(i)/kef;
// }

// void vychitanie(std::vector<double>& matrica, std::vector<double>& part_b, int n, int j)
// {
//     int i;
//     double kef;

//     for(i = 0; i < n; i++)
//     {
//         if(i == j)
//         {
//             continue;
//         }

//         kef = A(i, j);
//         //fprintf(stdout,"i=%d j=%d A(i,j)=%e\n",i,j,A(i,j));
        
//         for(int k = 0; k <n; k++)
//         {
//             A(i, k) = A(i,k) - A(j, k) * kef;
//         }
//         B(i) = B(i) - B(j) * kef;
//     }
// }

// double norma_matrica(std::vector<double>& matrica, int n)
// {
//     double sum,max;
//     max=0.;
//     for(int i=0;i<n;i++)
//     {
//         sum=0.;
//         for(int j=0;j<n;j++)
//         {
//             //fprintf(stdout,"i=%d j=%d A(i,j)=%e\n",i,j,A(i,j));
//             sum=sum+fabs(A(i,j));
//         }
        
//         if(sum>max)
//         {
//             max=sum;
//         }
//     }
//     return max;
// }

