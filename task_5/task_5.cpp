#include "task_5.hpp"


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


double f14(double x) { return sin(x); }
double f15(double x) { return cos(x); }
double f16(double x) { return tan(x); }
double f17(double x) { return log(x); }
double f18(double x) { return x * x - 4 * x + 3; }
double f19(double x) { return x * x * x - 3 * x * x + 2 * x; }
double f20(double x) { return x * sin(x); }
double f21(double x) { return exp(-x * x); }
double f22(double x) { return sqrt(x); }
double f23(double x) { return 1 / x; }
double f24(double x) { return exp(x) - 1; }
double f25(double x) { return (x != 0) ? sin(x) / x : 1; }  // Обработка деления на ноль
double f26(double x) { return x * cos(x); }
double f27(double x) { return x * x * x * x; }
double f28(double x) { return sqrt(1 - x * x); }
double f29(double x) { return log(1 + x); }
double f30(double x) { return x * x * exp(x); }


double g0(double x){ return -10. + x * 0;}
double g1(double x){ return -11.*2. * x;}
double g2(double x){ return -12.*3. * x*x;}
double g3(double x){ return -13.*4 * x*x*x*x;}
double h0(double x){ return 1 + x * 0;}
double h1(double x){ return x ;}
double h2(double x){ return x * x;}
double h3(double x){ return x * x * x;}
double h4(double x){ return x * x * x * x;}
double h5(double x){ return x * x * x * x * x;}
double h6(double x){ return x * x * x * x * x * x;}
double Runge(double x) {return 1 / (1 + 25 * x * x); }

double FunctionWithName::TrueRes(){
    double trueRes = 0;
    if (id < 4) {
        trueRes = pow(this->b, id + 1) - pow(this->a, id + 1);
    }
    //  f4 (f(x) = e^x)
    else if (id == 4) {
        trueRes = exp(b) - exp(a);
    }
    // f5 (f(x) = x * e^x)
    else if (id == 5) {
        //std::cout << "HERE!" << std::endl;
        trueRes = exp(b) * (b - 1) - exp(a) * (a - 1);
    }
    //  f6 (f(x) = cos(100x))
    else if (id == 6) {
        trueRes = (sin(100 * a) - sin(100 * b)) / 100.;
    }
    //  f7 (f(x) = e^(-10x))
    else if (id == 7) {
        trueRes = (exp(a) - exp(b)) / 1000.;
    }
    // f8 (f(x) = 1/sqrt(x))
    else if (id == 8) {
        trueRes = 2 * (sqrt(b) - sqrt(a));
    }
    
    else if (id > 9 && id < 14) {
        trueRes = id * (pow(this->a, this->id - 9) - pow(this->b, this->id - 9));
    }

    else if (id == 14) {  // f(x) = sin(x)
        trueRes = -cos(b) + cos(a);
    }
    else if (id == 15) {  // f(x) = cos(x)
        trueRes = sin(b) - sin(a);
    }
    else if (id == 16) {  // f(x) = tan(x)
        trueRes = log(cos(a)) - log(cos(b));  // log(cos(x)) — интеграл от tan(x)
    }
    else if (id == 17) {  // f(x) = log(x)
        trueRes = log(b) - log(a);
    }
    else if (id == 18) {  // f(x) = x^2 - 4x + 3
        trueRes = (pow(b, 3) / 3 - 2 * pow(b, 2) + 3 * b) - (pow(a, 3) / 3 - 2 * pow(a, 2) + 3 * a);
    }
    else if (id == 19) {  // f(x) = x^3 - 3x^2 + 2x
        trueRes = (pow(b, 4) / 4 - pow(b, 3) + pow(b, 2)) - (pow(a, 4) / 4 - pow(a, 3) + pow(a, 2));
    }
    else if (id == 20) {  // f(x) = x * sin(x)
        trueRes = -b * cos(b) + a * cos(a) + sin(b) - sin(a);
    }
    else if (id == 21) {  // f(x) = e^(-x^2)
        trueRes = (exp(-pow(b, 2)) - exp(-pow(a, 2)));
    }
    else if (id == 22) {  // f(x) = sqrt(x)
        trueRes = (2.0 / 3.0) * (pow(b, 3.0 / 2.0) - pow(a, 3.0 / 2.0));
    }
    else if (id == 23) {  // f(x) = 1/x
        trueRes = log(b) - log(a);
    }
    else if (id == 24) {  // f(x) = exp(x) - 1
        trueRes = exp(b) - exp(a);
    }
    else if (id == 25) {  // f(x) = sin(x)/x
        trueRes = (sin(b) / b) - (sin(a) / a);
    }
    else if (id == 26) {  // f(x) = x * cos(x)
        trueRes = sin(b) + b * cos(b) - (sin(a) + a * cos(a));
    }
    else if (id == 27) {  // f(x) = x^4
        trueRes = (pow(b, 5) / 5) - (pow(a, 5) / 5);
    }
    else if (id == 28) {  // f(x) = sqrt(1 - x^2)
        // Это определённый интеграл для арксинуса, что приводит к стандартному решению
        trueRes =  (b * sqrt(1 - b * b) + asin(b) - a * sqrt(1 - a * a) - asin(a));
    }
    else if (id == 29) {  // f(x) = log(1 + x)
        trueRes = log(1 + b) - log(1 + a);
    }
    else if (id == 30){  // f(x) = x^2 * exp(x)
        trueRes = (b * b - a * a) * exp(b) - 2 * (exp(b) - exp(a));
    }
    else if (id == 31){  // f(x) = 1/(1 + 25x^2)
        trueRes = (atan(5*b))/5 - (atan(5*a))/5;
    }
    
    else {
        trueRes = M_PI;  
    }

    return trueRes;
}


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

void writeResultsToFile(const std::string& filename, std::vector<FunctionWithName>& functions, int N) {
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
                << std::setw(20) << CompositeRectangleMethod(functions[i].a, functions[i].b, functions[i].func, N)
                << std::setw(20) << CompositeSimpsonQuadrature(functions[i].a, functions[i].b, functions[i].func, N)
                << std::setw(20) << CompositeGaussianQuadrature(functions[i].a, functions[i].b, functions[i].func, N)
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
             << std::setw(15) << std::fixed << std::setprecision(10) << N
             << std::setw(20) << std::fixed << std::setprecision(10) << error_R
             << std::setw(20) << std::fixed << std::setprecision(10) << error_S
             << std::setw(20) << std::fixed << std::setprecision(10) << error_G
             << std::endl;

        N++; 
    }

    fout.close();
    

}