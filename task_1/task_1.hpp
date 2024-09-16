#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

double f(double x);
double scalar(double *x, double*y, unsigned int N);
int    NodesGenerate(double *grid, int N); //генерация узлов
double Phi(int m, int k, int N);  // функции, вычисляющие и записывающие базисный вектор 
void   WritePhiTo(int m, int N, double* ph);
int    CoeffCalculate(int N, double* c_m, double (*f)(double), double* x_k, double* u_k, double* phi);
double FourierCompute(double* coefs, int N, double x);
double FullCompute(double x, int N, double* c_m, double (*f)(double), double* nodes, double * u_k, double* phi);



// else
//     {
//         WriteResult(netmemory, N, cNks, u, fp);
//         if ((sk = fopen("printAll.gpi", "w+")) == NULL)
//         {
//             printf("Не удалось открыть файл");
            
//             fclose(fp);
//             return 0;
//         }
//         WriteSkrypt(N, argv[2], sk);
//         fclose(fp);
//         fclose(sk);
//     }