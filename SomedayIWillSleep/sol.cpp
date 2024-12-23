#include "sol.hpp"

void phiCalculating(double *Phi, int m, int N)
{
	double h = 1. / ((double)N - 1.);

	for(int k = 0; k < N+1; k++)
	{
		Phi[k] = cos(M_PI * 0.5 * (m - 1.) * (2*k - 1) * h);
		// std::cout << Phi[k] << std::endl;
	}
}

void printVector(double *Vector, int N){
	for(int i = 0; i < N; i++){
		std::cout << Vector[i] <<"  ";
	}
    std::cout<<std::endl;
}

double scalarProduct(double *Phi, double *ValuesInNodes, int N)
{
	double res = 0.,
		   norm = 0.,
		   h = 1. / ((double)N - 1.);

	for(int i = 1; i < N; i++)
	{
		res += Phi[i] * ValuesInNodes[i] * h;
		norm += Phi[i] * Phi[i] * h;
	}

    // res += Phi[N] * ValuesInNodes[N] * h * 0.5;
	// norm += Phi[N] * Phi[N] * h * 0.5;
	res = res / (double)norm;
	return res;
}

double coefCalculating(double *Coef, double *Phi, double *ValuesInNodes, int N)
{
	double c = 0.;
	for(int m = 1; m < N; m++)
	{
		phiCalculating(Phi, m, N);
		c = scalarProduct(Phi, ValuesInNodes, N);
		Coef[m] = c;
	}	
	return 0;
}

double backFourier(double *Coef, double x, int N)
{
	double res = 0.;
	for(int m = 1; m < N; m++)
	{
		res += Coef[m] * cos(M_PI * (m - 1) * x);
	}
	return res;
}

void aproximateValuesCalculating(double *Coef, double *AproximateValues, double *Nodes, int N)
{
	for (int k = 0; k <= N+1; k++)
	{
		AproximateValues[k] = backFourier(Coef, Nodes[k], N);
	}
}

void copyCoefToFile(int N, double *Coef){
	std::ofstream outFile("coef.txt");
	for (int i = 1; i < N; i++) {
        outFile << Coef[i] << " " << std::endl;
	}
    outFile.close();
}

void copyPointsToFile(int N, double *Nodes, double *ValuesInNodes, double *AproximateValues){
	std::ofstream outFile("out.txt");
	for (int i = 1; i < N; i++) {
        outFile << Nodes[i] << " " << ValuesInNodes[i] << " " << AproximateValues[i] << std::endl;
	}
    outFile.close();
}

double normFunction(double (*f)(double), double *Coef, int N){
    double h = 1/10000.;
    double max = -1;
    double delta = 0;
    for (double x = 0.; x < 1.; x += h){
        delta = fabs(f(x) - backFourier(Coef, x, N));
        if (delta > max)
            max = delta;
    }
    return max;
}



double Lambda(int n, int N, double p){
    //double lam = p - 2. * ((double)(N-1) * (double)(N-1)) * (cos(M_PI * (n-1.) / (double)(N-1)) - 1.);
    double lam = p + 4. * ((double(N-1) * double(N-1))) * sin(M_PI * (n-1) / (2. * (double)(N-1))) * sin(M_PI * (n-1) / (2. * (double)(N-1)));
    
    return lam;
}


void MatrixFill(double p, double* M, int N){

    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            if (i == j){
                M[i * (N + 1) + j] = p + 2.0 * (double)((N-1) * (N-1));
            } 
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N + 1) + j] = -(double)((N-1) * (N-1));
            } 
            else{
                M[i * (N + 1) + j] = 0.0;
            }
        }
    }
    M[1 * (N + 1) + 1] = 1.* (double)((N-1.) * (N-1.)) + p;
    M[(N-1)*(N+1)+(N-1)] = 1.* (double)((N-1.) * (N-1.)) + p;

    M[0 * (N + 1) + 0] = 1.;
    M[(N)*(N+1)+(N)]   = 1.;
    M[0 * (N + 1) + 1] = -1.;
    M[(N)*(N+1)+(N-1)] = -1.;
    M[(N-1)*(N+1)+(N)] = 0.;

}

void printMatrix(double* matrix, int N) {
    for (int i = 0; i < N+1; ++i) { 
        for (int j = 0; j < N+1; ++j) { 
            std::cout << std::setw(10) << std::setprecision(4) << matrix[i * (N+1) + j] << " ";
        }
        std::cout << std::endl;
    } 
}


void searchCoef(double *coef, double *b, double *phi, double p, int N)
{
	double cn = 0., val = 1.;

	// Вычисляем базисные функции,
	// применяем формулу для поиска коэффициентов
	for(int n = 1; n < N; n++)
	{
		phiCalculating(phi, n, N);
	
		cn = scalarProduct(phi, b, N);

		//printf("coef_%d = %lf ", n, cn);

		val = Lambda(n, N, p);
        //std::cout   << "val = " << val << std::endl;

		if (val) coef[n] = cn / val;
	}
}

void searchSol(double *coef, double *x, int N)
{
	x[0] = 0.;

	double h = 1. / (double)(N-1);

	for(int k = 1; k < N; k++)
	{
		x[k] = 0.;
		
		//phiCompute(phi, m, N);

		for(int m = 1; m < N; m++)
		{
			x[k] += coef[m] * cos(M_PI * 0.5 * (m - 1.) * (2*k - 1) * h);
		}
	}
}

void MultiplyMatrixByVector(double* matrix, double* res, double* vec, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            res[i] += matrix[i*N+j] * vec[j];
        }
    }
}