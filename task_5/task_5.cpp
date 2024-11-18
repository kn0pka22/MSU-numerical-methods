#include "task_5.hpp"

// double f(double x){ 
//     return 1 + x * 0;
//     // return x * exp(x);
//     // return 2 * x;
//     // return 3 * x * x;
//     // return 4 * x * x * x;
//     // return exp(x);
//     // return cos(100*x);         // I[0, PI] = 0
//     // return exp(-1000 * x);     // I[0,1] = 10^(-3) 
//     // return 1 / sqrt(x);
//     // return 1 / sqrt(1 - x * x);// I[-1, 1] = PI
// }

double f0(double x){ return 1 + x * 0;}
double f1(double x){ return 2 * x;}
double f2(double x){ return 3 * x * x;}
double f3(double x){ return 4 * x * x * x;}
double f4(double x){ return exp(x);}
double f5(double x){ return x * exp(x);}
double f6(double x){ return cos(100*x);} // I[0, PI] = 0
double f7(double x){ return exp(-10 * x);} // I[0,1] = 10^(-3) 
double f8(double x){ return  1 / sqrt(x);}  //return (std::fabs(x) < 1e-15) ? 1e-10 : 1 / sqrt(x); 
double f9(double x){ return 1 / sqrt(1 - x * x); } // I[-1, 1] = PI return ((std::fabs(x)-1) < 1e-15) ? 1e-10 : 1 / sqrt(1 - x * x); 
double g0(double x){ return -10. + x * 0;}
double g1(double x){ return -11.*2. * x;}
double g2(double x){ return -12.*3. * x*x;}
double g3(double x){ return -13.*4 * x*x*x*x;}

double FunctionWithName::TrueRes(){
    double trueRes = 0;
    if (id < 4) {
        trueRes = pow(this->b, id + 1) - pow(this->a, id + 1);
    }
    else if (id == 4) {
        trueRes = exp(b) - exp(a);
    }
    else if (id == 5) {
        trueRes = exp(b) * (b - 1) - exp(a) * (a - 1);
    }
    else if (id == 6) {
        trueRes = (sin(100 * a) - sin(100 * b)) / 100.;
    }
    else if (id == 7) {
        trueRes = (exp(a) - exp(b)) / 1000.;
    }
    else if (id == 8) {
        trueRes = 2 * (sqrt(b) - sqrt(a));
    }
    else if (id > 9) {
        //std::cout<<"HERE!"<<std::endl;
        trueRes = id * (pow(this->a, this->id- 9) - pow(this->b, this->id- 9));
    }
    else {
        trueRes = M_PI;  
    }
    return trueRes;
}


// bool isNumber(const std::string& str) {
//     std::regex re("^[+-]?[0-9]*\\.?[0-9]+$");  // Регулярное выражение для числа
//     return std::regex_match(str, re);
// }

double RectangleMethod(double a, double b, const std::function<double(double)>& f){
    return (b - a) * f((a + b) / 2.);
}

double SimpsonMethod(double a, double b, const std::function<double(double)>& f){
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}


double GaussMethod(double a, double b, const std::function<double(double)>& f) {
    double x0 = (a + b) / 2;
    double xminus = x0 - (b - a) / 2. * sqrt(3. / 5.);
    double xplus = x0 + (b - a) / 2. * sqrt(3. / 5.);

    return ((b - a) / 18.) * (5. * f(xminus) + 8. * f(x0) + 5. * f(xplus));
}

double CompositeRectangleMethod(double a, double b,const std::function<double(double)>& f, int N){
    double step = (b - a) /(double)N;
    double ans = 0.;

    for(int i = 0; i < N; ++i){
        ans += RectangleMethod(a + i * step, a + (i + 1) * step, f);
    } 

    return ans;
}

double CompositeSimpsonQuadrature (double a, double b, const std::function<double(double)>& f, int N) {
    double step = (b - a) /(double)N;
    double ans = 0;

    for(int i = 0; i < N; ++i){
        ans += SimpsonMethod(a + i * step, a + (i + 1) * step, f);
    } 

    return ans;
}

double CompositeGaussianQuadrature(double a, double b, const std::function<double(double)>& f, int N) {
    double step = (b - a) /(double)N;
    double ans = 0;

    for(int i = 0; i < N; ++i) {
        ans += GaussMethod(a + i * step, a + (i + 1) * step, f);
    } 

    return ans;
}



void writeResultsToFile(const std::string& filename, const std::vector<FunctionWithName>& functions, std::vector<double>& WolframResults) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }


    outFile << std::left << std::setprecision(3) << std::fixed;

    outFile << std::setw(20) << "Segment" 
            << std::setw(30) << "Function" 
            << std::setw(20) << "Rectangle" 
            << std::setw(20) << "Simpson" 
            << std::setw(20) << "Gauss" 
            << std::setw(20) << "Wolfram" << std::endl; 
    
    outFile << std::string(130, '-') << std::endl;

    for (size_t i = 0; i < functions.size(); ++i) {
        outFile << std::setw(1) << "[ " 
                << std::setw(6) << functions[i].a 
                << std::setw(2) << ", " 
                << std::setw(6) << functions[i].b 
                << std::setw(1) << "]"
                << std::setw(3) << " "
                << std::setw(30) << functions[i].name
                << std::setw(20) << RectangleMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << SimpsonMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << GaussMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << WolframResults[i] << std::endl;
    }

    outFile.close();
}

void writeResultsToFile(const std::string& filename, std::vector<FunctionWithName>& functions) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }


    outFile << std::left << std::setprecision(3) << std::fixed;

    outFile << std::setw(20) << "Segment" 
            << std::setw(30) << "Function" 
            << std::setw(20) << "Rectangle" 
            << std::setw(20) << "Simpson" 
            << std::setw(20) << "Gauss" 
            << std::setw(20) << "True result" << std::endl; 
    
    outFile << std::string(130, '-') << std::endl;

    for (size_t i = 0; i < functions.size(); ++i) {
        outFile << std::setw(1) << "[ " 
                << std::setw(6) << functions[i].a 
                << std::setw(2) << ", " 
                << std::setw(6) << functions[i].b 
                << std::setw(1) << "]"
                << std::setw(3) << " "
                << std::setw(30) << functions[i].name
                << std::setw(20) << RectangleMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << SimpsonMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << GaussMethod(functions[i].a, functions[i].b, functions[i].func)
                << std::setw(20) << functions[i].TrueRes() << std::endl;
    }

    outFile.close();
}





void GenereteFileForPCalculation(int numTests, int numFunc, int N, const std::string& filename,std::vector<FunctionWithName>& functions) {
     double R, S, G, error_R, error_S, error_G;

    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }


    fout << std::left 
         << std::setw(15) << "N" 
         << std::setw(20) << "error_R" 
         << std::setw(20) << "error_S" 
         << std::setw(20) << "error_G" 
         << std::endl;


    for(int i = 0; i < numTests; ++i) {
        R = CompositeRectangleMethod(functions[numFunc].a, functions[numFunc].b, functions[numFunc].func, N);
        S = CompositeSimpsonQuadrature(functions[numFunc].a, functions[numFunc].b, functions[numFunc].func, N);
        G = CompositeGaussianQuadrature(functions[numFunc].a, functions[numFunc].b, functions[numFunc].func, N);

        error_R = (fabs(R - functions[numFunc].TrueRes()));
        error_S = (fabs(S - functions[numFunc].TrueRes()));
        error_G = (fabs(G - functions[numFunc].TrueRes()));

        
        fout << std::left 
             << std::setw(15) << std::fixed << std::setprecision(5) << N
             << std::setw(20) << std::fixed << std::setprecision(5) << error_R
             << std::setw(20) << std::fixed << std::setprecision(5) << error_S
             << std::setw(20) << std::fixed << std::setprecision(5) << error_G
             << std::endl;

        N *= 2; 
    }

    fout.close();
    

}