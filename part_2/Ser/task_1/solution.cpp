#include "solution.hpp"

double f(double x){
    // return cos(M_PI * 0. * x); 
    // return cos(M_PI * 1 * x);
    // return cos(M_PI * 2 * x);
    // return cos(M_PI * 3 * x);
    // return cos(M_PI * 4 * x);
    // return cos(M_PI * 5 * x);
    // return cos(M_PI * 0 * x) + cos(M_PI * 5 * x);
    // return 1 + 0. * x;   
	// return fabs(x)-1.;  
    // return x*x*x*exp(x);
    // return x*x;
    return 2./3. * x*x*x - x*x;
    // return 2./3. * pow(x, 12) - pow(x, 8);
    //return x*x * (1-x)*(1-x) * exp(-x);
}

// Задаём ВСЕ узлы
void fillingNodes(double *Nodes, int N){
    // x_0, ..., x_N
    for(int i = 0; i < N+1; i++){
		double h = 1. / (N - 1);
        Nodes[i] = -h/2 + i * h;
    }
}

// Задаём только ВНУТРЕННИЕ узлы
void fillingBasicNodes(double *BasicNodes, int N){
    double h = 1./(double)(N-1);
    // x_1, .., x_N-1 
    for(int i = 0; i < N-1; i++){
        BasicNodes[i] = h/2. + i * h;
    }
}

// Задаём значения во ВСЕХ узлах
void fillingValuesInNodes(double *Nodes, double *ValuesInNodes, int N){
    for(int i = 0; i < N+1; i++){
        ValuesInNodes[i] = f(Nodes[i]);
    }
}

// Задаём значения только во ВНУТРЕННИХ узлах
void fillingValuesInBasicNodes(double *BasicNodes, double *ValuesInBasicNodes, int N){
    for(int i = 0; i < N-1; i++){
        ValuesInBasicNodes[i] = f(BasicNodes[i]);
    }
}

// Заполняем матрицу B
void fillingMatrixB(double *B, double p, int N){
    for(int i = 0; i < N-1; i++){
        for(int j = 0; j < N-1; j++){
            if (i == j){
                B[i * (N - 1) + j] = p + 2.0 * (double)((N-1) * (N-1));
            } 
            else if ((i - j == 1 || i - j == -1)){
                B[i * (N - 1) + j] = -(double)((N-1) * (N-1));
            } 
            else{
                B[i * (N - 1) + j] = 0.0;
            }
        }
    }
    // k=1: (p+1/h^2) * y_1 - 1/h^2 * y_2
    B[0] = 1. * (double)((N-1) * (N-1)) + p;

    // - 1/h^2 * y_(N-2) + (p+1/h^2) * y_(N-1)
	B[(N-2)*(N-1)+(N-2)] = 1. * (double)((N-1) * (N-1)) + p;   
}

// Заполняем матрицу A
void fillingMatrixA(double *A, double p, double *pk, int N){
    for(int i = 0; i < N-1; i++){
        for(int j = 0; j < N-1; j++){
            if (i == j){
                A[i * (N - 1) + j] = p + pk[i] + 2.0 * (double)((N-1) * (N-1));
            } 
            else if ((i - j == 1 || i - j == -1)){
                A[i * (N - 1) + j] = -(double)((N-1) * (N-1));
            } 
            else{
                A[i * (N - 1) + j] = 0.0;
            }
        }
    }
    // k=1: (p+1/h^2) * y_1 - 1/h^2 * y_2
    A[0] = 1. * (double)((N-1) * (N-1)) + p + pk[0];

    // - 1/h^2 * y_(N-2) + (p+1/h^2) * y_(N-1)
	A[(N-2)*(N-1)+(N-2)] = 1. * (double)((N-1) * (N-1)) + p + pk[N-2];   
}

// Заполняем полную матрицу Фурье (N+1)*(N+1)
// void fillingFullMatrix(double p, double *Matrix, int N){
//     for(int i = 0; i < N+1; i++){
//         for (int j = 0; j < N+1; j++){
//             if (i == j){
//                 Matrix[i * (N + 1) + j] = p + 2.0 * (double)((N-1) * (N-1));
//             } 
//             else if ((i - j == 1 || i - j == -1)){
//                 Matrix[i * (N + 1) + j] = -(double)((N-1) * (N-1));
//             } 
//             else{
//                 Matrix[i * (N + 1) + j] = 0.0;
//             }
//         }
//     }
//     // y_0 = y_1
//     Matrix[0 * (N+1) + 0]= 1.;
//     Matrix[0 * (N+1) + 1]= -1.;

//     // k=1: 0 * y_0 + (p+1/h^2) * y_1 - 1/h^2 * y_2
//     Matrix[1 * (N+1) + 0]= 0.;
//     Matrix[1 * (N+1) + 1]= p + 1. * (double)((N-1) * (N-1));

//     // y_(N-1) = y_N
//     Matrix[N * (N+1) + N]= -1.;
//     Matrix[N * (N+1) + N-1]= 1.;

//     // k = N-1: - 1/h^2 * y_(N-2) + (p+1/h^2) * y_(N-1) + 0 * y_N
//     Matrix[(N-1)* (N+1) + N-1]= p + 1. * (double)((N-1) * (N-1));
//     Matrix[(N-1) * (N+1) + N]= 0.;
// }

// Умножение вектора на матрицу Фурье с постоянной добавкой p на диагонали (Умножение на матрицу B)
void MultiplicationByB(double *vector, double *res, int N, double p){
    for(int i = 1; i < N; i++){
        if (i == 1){
            res[i-1] = (p + double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i];
        }
        else if (i == N-1){
            res[i-1] = (p + double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i-2];
        }
        else{
            res[i-1] = (p + 2. *double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i-2] - double((N-1)*(N-1)) * vector[i];
        } 
    }
}

// Умножение вектора на матрицу Фурье с постоянной добавкой p и добавками pk на диагонали (Умножение на матрицу A)
void MultiplicationByA(double *vector, double *res, int N, double p, double *pk){


    //на вход получилиNn-1
    for (int i = 1; i < N; i++){
        if (i == 1){
            res[i-1] = (pk[i-1] + p + double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i];
        }
        else if (i == N-1){
            res[i-1] = (pk[i-1] + p + double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i-2];
        }
        else{
            res[i-1] = (pk[i-1] + p + 2. *double((N-1)*(N-1))) * vector[i-1]  - double((N-1)*(N-1)) * vector[i-2] - double((N-1)*(N-1)) * vector[i];
        }
    }
}

// Вычисляем собственные значения матрицы Фурье с постоянной добавкой p на диагонали
void Lambda(double *lambda, double p, int N){
    for(int m = 1; m < N; m++){
        lambda[m-1] = p + 4. * ((double(N-1) * double(N-1))) * sin(M_PI * (m-1) / (2. * (double)(N-1))) * sin(M_PI * (m-1) / (2. * (double)(N-1)));;
    }
}

// Вычисляем m-ый собственный вектор (функцию) 
void phiCalculating(double *Phi, int m, int N){
	double h = 1. / ((double)N - 1.);

	for(int k = 1; k < N; k++){
		Phi[k-1] = cos(M_PI * 0.5 * (m - 1.) * (2*k - 1) * h);
	}
}

/// Скалярное проивезедение (Phi, Vector)_h -> даёт соответствующий коэффициент (res) в разложении Vector по базису из собственных векторов
double scalarProduct(double *Phi, double *b, int N){
	double res = 0.,
		   norm = 0.,
		   h = 1. / ((double)N - 1.);

	for(int k = 1; k < N; k++){
		res += Phi[k-1] * b[k-1] * h; // Скалярное произведение в числителе
		norm += Phi[k-1] * Phi[k-1] * h; // Скалярное произведение в знаменателе
	}
	res = res / (double)norm;
	return res;
}


// Координаты вектора 
bool searchCoef(double *coef, double *b, double *phi, double p, int N, double *lambda){ // ?Bx = b?
	double dm = 0.;
    double eigenValue = 1.;
    Lambda(lambda, p, N);
    double dm1 = 0.;

	// Вычисляем базисные функции, применяем формулу для поиска коэффициентов
	for(int m = 1; m < N; m++){
		phiCalculating(phi, m, N);
	
		dm = scalarProduct(phi, b, N);
        if (m==1) dm1 = dm;

		eigenValue = lambda[m-1];
        //std::cout<<"m = "<<m<<" p = "<<p<<" dm = "<<dm<<" fabs(dm) = "<<fabs(dm)<<" eigenValue = "<<eigenValue<<std::endl;
        if (fabs(p)<1e-15){
            if (fabs(dm1)<1e-15){
                if (fabs(eigenValue)>1e-15) coef[m-1] = dm / eigenValue; 
                else coef[m-1] = 0.;
            }
            else{
                std::cout<<"Invalid requirement!  \\_/" << std::endl;
                return false;
            }            
        }
        else{
            if (eigenValue) coef[m-1] = dm / eigenValue; 
            else coef[m-1] = 0.;

        } 
        
        // if ((((fabs(p)<1e-15) && (m==1) && (fabs(dm)<1e-15))) || (fabs(p)>1e-15)){
        // //if (val){ 
        //     coef[m-1] = dm / eigenValue; 
        // }
		// else{ 
        //     std::cout<<"Invalid requirement!  \\_/" << std::endl; 
        //     return false;
        // } 
	}
    return true;
}

void searchSol(double *coef, double *x, int N){
	x[0] = 0.;

	double h = 1. / (double)(N-1);

	for(int k = 1; k < N; k++){
		x[k-1] = 0.;

		for(int m = 1; m < N; m++){
			x[k-1] += coef[m-1] * cos(M_PI * 0.5 * (m - 1.) * (2*k - 1) * h);
		}
	}
}


double BSolver( double* x, double* b, 
                double tau, int mIter, double* Ax, 
                double* memory1, double p, int N, double *pk,
                double *Coef, double *Phi, double *lambda){

    double mean = 0.;
  
    for (int k = 0; k < N-1; k++){
        x[k] = 0;
    }

    for (int iteration = 0; iteration < mIter; iteration++){
        
        MultiplicationByA(x, Ax, N, p, pk);  // Ax

        for (int k = 0; k < N-1; k++){ // b - Ax
            Ax[k] = b[k] - Ax[k];
        }

        for(int i = 0; i < N-1; i++){
            mean += (pk[i]+p);
        }
        
        mean /= ((double)(N-1));
    
        searchCoef(Coef, Ax, Phi, mean, N, lambda); // Решаем Фурье - находим y^(iteration+1)
        searchSol(Coef, memory1, N); // Решаем Фурье - находим y^(iteration+1) (B(y^(iteration+1)) = b - A(x^iteration)), memory1 = y^(iteration+1)

        for (int k = 0; k < N-1; k++){
            x[k] += tau * memory1[k]; // x^(iteration+1) = x^iteration + tau * y^(iteration+1)
        }  
    }
    return ErNormInf(b, x, memory1, N, pk, p);
}

// Норма вектора невязки: нашли x - приближенное решение уравнения Ax = b. Насколько действительно Ax отличается от b?
// double ErNormInf(double* b, double* x, double* memory, int N, double *pk, double p){
//     double maxDifference = 0.;
//     MultiplicationByA(x, memory, N, p, pk);
//     for (int k = 0; k < N-1; k++){
//         if (fabs((b[k] - memory[k])) > maxDifference){
//             maxDifference = fabs((b[k] - memory[k]));
//         }
//     }
//     return maxDifference;
// }


double ErNormInf(double* b, double* x, double* memory, int N, double *pk, double p){
    double residual_squared = 0.0;
    double h = 1.0 / (double)(N - 1);
    
    MultiplicationByA(x, memory, N, p, pk);
    
    for (int k = 0; k < N-1; k++){
        double diff = b[k] - memory[k];
        residual_squared += diff * diff;
    }
    
    return sqrt(residual_squared * h);
}

//==============================================================================================================================================//
// Печать матрицы в консоль
void printMatrix(double *Matrix, int N){
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			std::cout << std::setw(10) << std::setprecision(4) << Matrix[i * N + j] << " ";
		}
		std::cout << std::endl;
	}
}

// Печать вектора в консоль
void printVector(double *Vector, int N){
	for(int k = 0; k < N; k++){
		std::cout << std::setw(10)<< std::setprecision(15) << Vector[k] << " ";
	}
	std::cout << std::endl;
}

// Умножение вектора на собственное число (использовалось для проверки того, что собственные функции и собственные числа действительно собственные)
// void MultbyNum(double *vector, double *res, double *lambda, int N, int m){
//     for (int i = 0; i < N-1; i++){
//         res[i] = lambda[m-1] * vector[i];
//     }
// }
