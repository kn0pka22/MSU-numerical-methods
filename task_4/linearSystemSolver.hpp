#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <vector>
#define e(i, j, n) ((i - 1) * (n) + (j - 1))
#define EPS 1e-15


#define A(i, j) matrica[i * n + j]
#define B(i) part_b[i]



int solve(std::vector<double>& M, std::vector<double>& b, std::vector<double>& x, std::vector<int>& memory);
void printMatrix(const std::vector<double>& matrix);
void printVector(const std::vector<double>& vec);


// int Michael_Jordan(int n, std::vector<double>& matrica, std::vector<double>& part_b,std::vector<double>& result); //
// int lead(std::vector<double>& matrica, int n, int j, double z); //
// void delenie(std::vector<double>& matrica, std::vector<double>& part_b, int n, int i); //
// void vychitanie(std::vector<double>& matrica, std::vector<double>& part_b, int n, int i);//
// double norma_matrica(std::vector<double>& matrica, int n);//