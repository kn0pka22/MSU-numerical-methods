#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <vector>
#define e(i, j, n) ((i - 1) * (n) + (j - 1))



int solve(std::vector<double>& M, std::vector<double>& b, std::vector<double>& x, std::vector<int>& memory);
void printMatrix(const std::vector<double>& matrix);
void printVector(const std::vector<double>& vec);