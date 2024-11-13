
#include "solution.hpp"

double f(double x, double y){
    return cos(M_PI * x / 2.) * cos(M_PI * y / 2.);
    //return 1 + 0. * x + 0. * y;
    //return (2./3. * x*x*x - x*x) * (2./3. * y*y*y - y*y);
    //return 2./3. * pow(x, 12) - pow(x, 8);
    //return x*x * (1-x)*(1-x) * exp(-x); 
}

int main(){

    int N = 5;

    double *Nodes = new double[N+1];
    double *ValuesInNodes = new double[(N+1)*(N+1)];
    double *Phi = new double[N+1];
    double *Memory = new double[N+1];
    double *D = new double[(N+1)*(N+1)];
    double *C = new double[(N+1)*(N+1)];
    double *AproximateValues = new double[(N+1)*(N+1)];

    fillingNodes(Nodes, N);
    fillingValuesInNodes(Nodes, ValuesInNodes, N);
    //PrintMatrix(ValuesInNodes, N, "VIN");
    
    fillingDMatrix(D, Phi, ValuesInNodes, N);
    // PrintMatrix(D, N, "D");

    fillingCMatrix(C, Phi, Memory, D, N);
    // PrintMatrix(C, N, "C");
    // //printNodes(Nodes, N);

    coefCalculating(C, Phi, ValuesInNodes, N);


    double h = 1. / (N-1);
    h/=3;
    double xi = 1. / (double)(N-1);
    double yi = 1. / (double)(N-1);
    xi/=(-2);
    yi/=(-2);


    bool flag = 0;
    if (flag){
        for (int i=0;i<(N+1)/3.;++i){
            for (int j=0;j<(N+1)/3.;++j){
                std::cout<<" backFourier( "<<xi<<", "<<yi <<" ) = "<<backFourier(C,xi,yi,N)<<" and f( "<<xi<<", "<<yi <<" ) = "<<f(xi,yi)<<std::endl;  
                xi += h;
                yi += h;
                std::cout<<" backFourier( "<<xi<<", "<<yi <<" ) = "<<backFourier(C,xi,yi,N)<<" and f( "<<xi<<", "<<yi <<" ) = "<<f(xi,yi)<<std::endl;  
                xi += h;
                yi += h;
                std::cout<<" backFourier( "<<xi<<", "<<yi <<" ) = "<<backFourier(C,xi,yi,N)<<" and f( "<<xi<<", "<<yi <<" ) = "<<f(xi,yi)<<std::endl;  
                xi += h;
                yi += h;
            }
        }
        xi += h;
        yi += h;
        std::cout<<" backFourier( "<<xi<<", "<<yi <<" ) = "<<backFourier(C,xi,yi,N)<<" and f( "<<xi<<", "<<yi <<" ) = "<<f(xi,yi)<<std::endl;  
        xi += h;
        yi += h;
        std::cout<<" backFourier( "<<xi<<", "<<yi <<" ) = "<<backFourier(C,xi,yi,N)<<" and f( "<<xi<<", "<<yi <<" ) = "<<f(xi,yi)<<std::endl;  
    }
    else{
        std::ofstream outFile("out.txt");
            if (outFile.is_open()) {
                double h = 1. / (N-1);
                h/=3;
                double xi = 1. / (double)(N-1);
                double yi = 1. / (double)(N-1);
                xi/=(-2);
                yi/=(-2);

                outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
                <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
                <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
                <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;
                for (int i = 1; i < N-1; ++i){ 
                    for (int j = 1; j < N; ++j){ 
                        outFile << std::setprecision(15) << std::fixed \
                        << std::setw(20) << xi << " " \
                        << std::setw(20) << yi << " " \
                        << std::setw(20) << f(xi,yi) << " " \
                        << std::setw(20) << backFourier(C,xi,yi,N) << std::endl;

                        xi += h;   
                        yi += h;
                        outFile<< std::setprecision(15) << std::fixed \
                        << std::setw(20) << xi << " " \
                        << std::setw(20) << yi << " " \
                        << std::setw(20) << f(xi,yi) << " " \
                        << std::setw(20) << backFourier(C,xi,yi,N) << std::endl;

                        xi += h;   
                        yi += h;
                        outFile<< std::setprecision(15) << std::fixed \
                        << std::setw(20) << xi << " " \
                        << std::setw(20) << yi << " " \
                        << std::setw(20) << f(xi,yi) << " " \
                        << std::setw(20) << backFourier(C,xi,yi,N) << std::endl;
                    }
                }
                xi += h;   
                yi += h;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << backFourier(C,xi,yi,N) << std::endl;

                xi += h;   
                yi += h;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << backFourier(C,xi,yi,N) << std::endl;
            }
            else {
                std::cerr << "Error opening file" << std::endl;
            }
        }

    

    //WriteToConsole(N, Nodes, ValuesInNodes, C,D, Memory, Phi);

    //aproximateValuesCalculating(C, AproximateValues, Nodes, N);

    //copyCoefToFile(N, C);

    //copyPointsToFile(N, Nodes, ValuesInNodes, AproximateValues);

    //p(N);

    delete[] Nodes;
    delete[] ValuesInNodes;
    delete[] Phi;
    delete[] D;
    delete[] C;
    delete[] AproximateValues;
    return 1;
}
