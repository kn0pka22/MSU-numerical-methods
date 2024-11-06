#include "task_4.hpp"

double f(double x){
    //return((exp(x)-exp(1.0))*sin(x));
    return fabs(x);
    //return (x*x+sin(x))*cos(3*x);
    //return 1./(1.+25.*x*x); //Runge's function
    //return std::fabs(x);
}

int GenerateEquidistantNodes(double a, double b, std::vector<double>& nodes){
    int N = nodes.size();
    if (N < 2 || a >= b) {
        return -1; 
    }

    double delta = (b - a) / double(N - 1);
    double now = a;

    for (int i = 0; i < N; ++i) {
        nodes[i] = now;          
        //values[i] = f(now);    
        now += delta;         
    }
    return 0; 
}

int GenerateChebyshevNodes(double a, double b, std::vector<double>& nodes){
    int N = nodes.size();

    if (N < 1 || a >= b) {
        return -1; 
    }

    double delta1 = (b + a) / 2;
    double delta2 = (b - a) / 2; //[a,b] <-> [-1,1]

    for (int i = 0; i < N; ++i) {
        nodes[N - 1 - i] = delta1 + delta2 * cos((2 * i + 1) * M_PI / (2 * N)); 
        //values[N - 1 - i] = f(nodes[N - 1 - i]); 
    }
    return 0;     
}

void ExtendedNodes(std::vector<double>& ExNodes, std::vector<double>& Nodes){
    double h;
    int n = Nodes.size();
    //int nn = 3*n; 
    for (int i=0;i<n;i++){
        ExNodes[3*i]=Nodes[i];
    }
    for(int i=0;i<n-1;i++){
        h = (Nodes[i+1]-Nodes[i])/3;
        ExNodes[3*i+1]=ExNodes[3*i]+h;
        ExNodes[3*i+2]=ExNodes[3*i+1]+h;
    }
}

void ExtendedValues(std::vector<double>& ExValues, std::vector<double>& ExNodes){
    int n = ExNodes.size();
    for(int i=0;i<n;i++){
        ExValues[i]=f(ExNodes[i]);
    }
}

void FillingValues( std::vector<double>& nodes,  std::vector<double>& values, double (*f)(double), int N){

    for (int i=0;i<N+1;++i){
        values[i] = f(nodes[i]);
    }
}


void ExtendedF(std::vector<double>& ExF, const std::vector<double>& result, const std::vector<double>& ExNodes, int MM) {
    for (int i = 0; i < 3 * MM - 2; i++) {
        double ans = result[1]; 
        for (int j = 1; j < MM - 1; j++) {
            ans += result[j+1] * std::pow(ExNodes[i], j); 
        }
        ExF[i] = ans; 
    }
}


int MatrixFill(std::vector<double>& matrix, const std::vector<double>& nodes){
    int n = sqrt(matrix.size());
    for (int i = 0; i < n; ++i) { 
        for (int j = 1; j < n; ++j) { 
            matrix[i*n+j] = pow(nodes[i],(j-1));
        }
        matrix[i * n] = (i % 2 == 0) ? 1.0 : -1.0;
    }
    return 0;
}

double delta(std::vector<double>& ExValues, std::vector<double>& ExF, int MM){
    double max=0.;
    double del;
    for(int i=0;i<3*MM-2;++i){
        del = std::fabs(ExValues[i]-ExF[i]);
        if(del>max){
            max=del;
        }
    }
    return max;
}

void CreateSigma(std::vector<double>& sigma, std::vector<double>& nodes, int MM, int N){
    int k;
    k=N/MM;
    for(int i=0;i<MM;++i){
        sigma[i]=nodes[i*k];
    }
}

double MaxDeviation(const std::vector<double>& nodes, std::vector<double>& sigma, const std::vector<double>& res, const std::vector<double>& valuesAll, int MM, int N){
    int change = -1, b;
    double h = res[0]; 
    double max = std::fabs(h);

    //Search for max value
    for (int i = 0; i < MM - 1; i++) {
        double d = std::fabs(res[i]);
        if (d > max) {
            max = d;
            change = i;
        }
    }

    // check
    if ( std::fabs(h)+ 1e-4 > max) {
        return 0; 
    }
    else{
        if (nodes[change] < sigma[0]) {
            double a = FindApprox(res, sigma[0], MM);
            double m = FindApprox(res, nodes[change], MM);
            double c = f(sigma[0]);
            if (((a - c < 0) && (m - valuesAll[change] < 0) && (m - valuesAll[change] > 0))) {
                sigma[0] = nodes[change];
            } 
            else{
                for (int i = 0; i < MM - 1; i++) {
                    sigma[MM - 1 - i] = sigma[MM - 2 - i];
                }
                sigma[0] = nodes[change];
            }
        } else if (nodes[change] > sigma[MM - 1]) {
            double a = FindApprox(res, sigma[MM - 1], MM);
            double m = FindApprox(res, sigma[change], MM);
            double c = f(sigma[MM - 1]);
            if (((a - c < 0) && (m - valuesAll[change] < 0) && (m - valuesAll[change] > 0))) {
                sigma[MM - 1] = nodes[change];
            } 
            else {
                for (int i = 0; i < MM - 1; i++) {
                    sigma[i] = sigma[i + 1];
                }
                sigma[MM - 1] = nodes[change];
            }
            return -1;
        } 
        else {
            b = 0;
            while (sigma[b] < nodes[change]) {
                b = b + 1;
            }
            double a = FindApprox(res, sigma[b - 1], MM);
            double m = FindApprox(res, nodes[change], MM);
            double c = f(sigma[b - 1]);
            double q = f(nodes[change]);
            if (((a - c < 0) && (m - q < 0) && (m - q > 0))) {
                sigma[b - 1] = nodes[change];
            } else {
                sigma[b] = nodes[change];
            }
            return -1;
        }
    }
}


double FindApprox(const std::vector<double>& result, double x, int z) {
    double ans = 0.0;

    for (int i = 0; i < z; i++) {
        int power = z - 1 - i; 

        double st = (power > 0) ? 1.0 : 0.0;

        // Если степень больше 0, вычисляем st = x^power
        if (power > 0) {
            for (int j = 0; j < power; j++) {
                st *= x; // Умножаем st на x
                if (std::fabs(st) < 1e-15) {
                    st = 0.0; 
                    break; 
                }
            }
        }

        ans += st * result[i];
    }

    return ans; // Возвращаем результат
}
