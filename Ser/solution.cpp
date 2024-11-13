#include "solution.hpp"
#include <iomanip>

void fillingNodes(double *Nodes, int N){
    double h = 1. / (N-1);
    for (int i = 0; i < N+1; i++){
        Nodes[i] = -h/2 + i*h;
    }
}
void printNodes(double *Nodes, int N){
	std::ofstream file("Nodes.txt");
	for (int i = 0; i < N+1; i++){
		file << std::setprecision(15) << Nodes[i] << std::endl;
	}
	file.close();
}
// void printValueInNodes(double *ValuesInNodes, int N){
// 	std::ofstream file("ValuesInNodes.txt");
// 	for (int i = 0; i < N+1; i++){
// 		file << std::setprecision(15) << ValuesInNodes[i] << std::endl;
// 	}
// 	file.close();
// }

// void printMatrix(double *matrix, int N){
//     for (int i = 0; i < N+1; ++i) { // Перебор строк
//         for (int j = 0; j < N+1; ++j) { // Перебор столбцов
//             std::cout << std::fixed << std::setprecision(10) << matrix[i * (N+1) + j] << " "; // Вывод элемента
//         }
//         std::cout << std::endl; // Переход на новую строку после вывода строки матрицы
//     }
// }

void fillingValuesInNodes(double *Nodes, double *ValuesInNodes, int N){
    for (int i = 0; i < N+1; i++){
		for (int j = 0; j < N+1; j++){
			ValuesInNodes[i * (N + 1) + j] = f(Nodes[i], Nodes[j]);
		}
    }
}

void phiCalculating(double *Phi, int m, int N)
{
	double h = 1. / ((double)N - 1.);

	for (int k = 0; k < N+1; k++)
	{
		Phi[k] = cos(M_PI * 0.5 * (m - 1.) * (2*k - 1) * h);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void printMatrix(double *Matrix, int N){
	std::ofstream file("Nodes.txt");
	for (int i = 0; i < N+1; i++){
		for (int j = 0; j < N+1; j++){
			
		}
		file << std::setprecision(15) << Nodes[i] << std::endl;
	}
	file.close();
}*/

void PrintMatrix(const double* matrix, int N, const std::string& name) {
    std::cout << "==================================================================" << name << "==================================================================" << std::endl;

    for (int k = 0; k < N * N; k++) {
        if (k % N == 0 && k != 0) std::cout << std::endl;
        std::cout << std::fixed << std::setprecision(15) << std::setw(20) << *(matrix + k) << "    "; 
    }
    std::cout << std::endl << std::endl;
}


double scalarProduct(double *Array1, double *Array2, int N) //Array1 = Phi, Array2 = ValuesInNodes
{
	double res = 0., norm = 0., h = 1. / ((double)N - 1.);
	for (int i = 1; i < N; i++)
	{
		res += Array1[i] * Array2[i] * h;
		norm += Array1[i] * Array1[i] * h;
	}
	res = res / (double)norm;
	return res;
}

void coefCalculating(double *Coef, double *Phi, double *Array1d, int N)
{
	double c = 0.;
	for (int m = 1; m < N; m++)
	{
		phiCalculating(Phi, m, N);
		c = scalarProduct(Phi, Array1d, N);
		Coef[m] = c;
	}	
}

void fillingDMatrix(double* D, double* Phi, double* ValuesInNodes, int N){
    for (int i = 1; i < N; i++)
	{ 
        coefCalculating(D + i * (N+1), Phi, ValuesInNodes + i * (N+1), N);
    }
	//PrintMatrix(D,N,"D");
}

void fillingCMatrix(double *C, double *Phi, double *Memory, double *D, int N){
	for (int j = 1; j < N; j++)
	{
        for (int i = 0; i < N+1; i++){
            Memory[i] = D[i *(N+1) + j];
			//std::cout<<"Memory[i] = "<<Memory[i]<<std::endl;
        }
        coefCalculating(C + j * (N+1), Phi, Memory, N);  
    }

}
double backFourier(double *C, double x, double y, int N)
{
	double res = 0.;
	for (int i = 1; i < N; i++)
	{
		for (int j = 1; j < N; j++){
			res += C[i*(N+1) + j] * cos(M_PI * (i - 1) * x) * cos(M_PI * (j - 1) * y);
		}
	}
	return res;
}

double WriteToConsole(int N, double* Nodes,  double* ValuesInNodes, double* C, double* D, double* Memory, double* Phi){
    // double h = 1/(N-0.5);
    double xi = 0;
    double yi = 0;
    double deltax = 0;
    double deltay = 0;
    std::cout<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
    <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
    <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
    <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl; 

    clock_t start=clock();
    fillingNodes(Nodes, N);
	printNodes(Nodes, N);

    fillingValuesInNodes(Nodes, ValuesInNodes, N);
	//PrintMatrix(ValuesInNodes, N, "VIN");

    fillingDMatrix(D, Phi, ValuesInNodes, N);
    //PrintMatrix(D, N, "D");

    fillingCMatrix(C, Phi, Memory, D, N);
    //PrintMatrix(C, N, "C");
    clock_t end=clock();
    double duration =(double)(end-start)/CLOCKS_PER_SEC;
    

    for (int i = 0; i < N+1; ++i){ 
        for (int j = 0; j < N+1; ++j){ 
            xi = Nodes[i];
            deltax = Nodes[i+1] - xi;
            deltax /= 3.;    
            yi = Nodes[j];
            deltay = Nodes[j+1] - yi;
            deltay /= 3.;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << backFourier(C, xi, yi, N) << std::endl;

            xi += deltax;   
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << backFourier(C, xi, yi, N) << std::endl; 

            xi += deltax;   
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << backFourier(C, xi, yi, N) << std::endl; 
        }
    }
    return duration;
}

/*void aproximateValuesCalculating(double* Coef, double *AproximateValues, double *Nodes, int N){
	for (int k = 0; k <= N+1; k++)
	{
		AproximateValues[k] = backFourier(Coef, Nodes[k], N);
	}
}*/

/*void copyCoefToFile(int N, double *Coef){
	std::ofstream outFile("coef.txt");
	for (int i = 1; i < N; ++i) {
        outFile << Coef[i] << " " << std::endl;
	}
    outFile.close();
}*/

/*void copyPointsToFile(int N, double *Nodes, double *ValuesInNodes, double *AproximateValues){
	std::ofstream outFile("out.txt");
	for (int i = 1; i < N; ++i) {
        outFile << Nodes[i] << " " << ValuesInNodes[i] << " " << AproximateValues[i] << std::endl;
	}
    outFile.close();
}*/

/*double normFunction(double (*f)(double), double *Coef, int N){
    double h = 1/10000.;
    double max = -1;
    double delta = 0;
    for (double x = 0.; x < 1.; x += h){
        delta = fabs(f(x) - backFourier(Coef, x, N));
        if (delta > max)
            max = delta;
    }
    return max;
}*/

