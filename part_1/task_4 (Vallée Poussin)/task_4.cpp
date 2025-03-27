#include "task_4.hpp"

double f(double x){
    //return((exp(x)-exp(1.0))*sin(x));
    return fabs(x);
    //return exp(x);
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

    for (int i=0;i<N;++i){
        values[i] = f(nodes[i]);
    }
}


void ExtendedF(std::vector<double>& ExF, const std::vector<double>& result, const std::vector<double>& ExNodes, int MM) {
    for (int i = 0; i < 3 * MM - 2; i++) {
        // double ans = result[1]; 
        // for (int j = 1; j < MM - 1; j++) {
        //     ans += result[j+1] * std::pow(ExNodes[i], j); 
        // }
        // ExF[i] = ans; 
        ExF[i] = CalcPolynom(result, ExNodes[i], MM);

    }
    std::ofstream outFile("data.txt");
    if (outFile.is_open()) {
        
    
        double tmp1,tmp2,err;

        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x)"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"Pn"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"err"<<std::setw(9)<<" "<<std::endl;
        
        int N = ExNodes.size();
        
        for (int i = 0; i < N; ++i){ 
            tmp1 = f(ExNodes[i]);
            tmp2 =  CalcPolynom(result,ExNodes[i],MM);
            err = std::fabs(tmp1-tmp2);
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << ExNodes[i]   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

        }
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
    //std::cout<<"FIRST CREATED SIGMA: "; printVector(sigma);
}

bool MaxDeviation(const std::vector<double>& nodes, std::vector<double>& sigma, const std::vector<double>& coeffs, std::vector<double>& values, const std::vector<double>& valuesAll, int MM, int N){
    double h = coeffs[0]; 
    //std::cout<<"h   = "<<h<<std::endl;

    double maxx = std::fabs(h);
    double res = 0.;
    double res_old, res_new;
    int k=0;
    
    //std::cout<<"Nodes now: \n"; printVector(nodes);
    
    //We are looking for the maximum deviation in the nodes
    for(int i=0;i<N;++i){  //phi
        //std::cout<<"CalcPolynom(coeffs, nodes[i], MM) = "<<CalcPolynom(coeffs, nodes[i], MM)<<std::endl;
        res = std::fabs(CalcPolynom(coeffs, nodes[i], MM) - valuesAll[i]);
        if (res > maxx){
            maxx = std::fabs(res);
            k=i;
           
            //std::cout<<"k=i = "<<k<<std::endl;
        }
    }
  
    //std::cout<<"max = "<<maxx<<std::endl;

    //exit and win!
    if (std::fabs(h) + 1e-6 > maxx){   
        //std::cout<<"fabs(h) + 1e-4 > maxx"<<std::endl;
		return 0;
	}
    else{
        //std::cout<<"in else"<<std::endl;
		if(nodes[k] < sigma[0]){
            //std::cout<<"case LEFT"<<std::endl;

			res_old = CalcPolynom(coeffs, sigma[0], MM);
			res_new = CalcPolynom(coeffs, nodes[k], MM);
            //std::cout<<"k = "<<k<<" and RES_NEW = "<<res_new<<std::endl;
            //printVector(sigma);
			
            if(((res_old - values[0] < 0) && (res_new - valuesAll[k] < 0)) || ((res_old - values[0] > 0) && (res_new - valuesAll[k] > 0))){
				sigma[0] = nodes[k];
				values[0] = valuesAll[k]; 
			}
			else{
				// We go in reverse order, because otherwise the values ​​are overwritten
                //printVector(sigma);
				for(int i = 0; i < MM - 1; ++i){
					sigma[MM - i - 1] = sigma[MM - 2 - i];
					values[MM - i - 1] = values[MM - 2 - i]; 
				}


				sigma[0] = nodes[k];
				values[0] = valuesAll[k]; 

			}
            //printVector(sigma);
			return 1;
		}
        else if(nodes[k] > sigma[MM - 1]){
            //std::cout<<"case RIGHT"<<std::endl;

			res_old = CalcPolynom(coeffs, sigma[MM - 1], MM);
			res_new = CalcPolynom(coeffs, nodes[k], MM);

            //printVector(sigma);
           

			if(((res_old - values[MM-1] < 0) && (res_new - valuesAll[k] < 0) ) || ( (res_old - values[MM-1] > 0) && (res_new - valuesAll[k] > 0))){
				sigma[MM - 1] = nodes[k];
				values[MM - 1] = valuesAll[k]; 
                //std::cout<<"HHHHHEEEEEEEERRRRRRRRREEEEEEEE\n";
			}
			else{
				for(int i = 0; i < MM - 1; ++i){
					sigma[i] = sigma[i + 1];
					values[i] = values[i + 1]; 
				}
				sigma[MM - 1] = nodes[k];
				values[MM - 1] = valuesAll[k]; 
			}
            //printVector(sigma);

			return 1;
		}
        else{
            //std::cout<<"case MIDDLE"<<std::endl;
            //printVector(values);
			int b = 0;
			while(sigma[b] < nodes[k]) b++;
            //std::cout<<"b = "<<b<<" and MM-1 = "<<MM-1<<std::endl;

			res_old = CalcPolynom(coeffs, sigma[b - 1], MM);
			res_new = CalcPolynom(coeffs, nodes[k], MM);

            // std::cout<<"res_old = "<<res_old<<std::endl;
            // std::cout<<"res_new = "<<res_new<<std::endl;
            // std::cout<<"values[b-1] = "<<values[b-1]<<std::endl;
            // std::cout<<"valuesAll[k] = "<<valuesAll[k]<<std::endl;

            //printVector(sigma);
            
			if(((res_old - values[b-1] < 0) && (res_new - valuesAll[k] < 0) ) || ( (res_old - values[b-1] > 0) && (res_new - valuesAll[k] > 0))){
				sigma[b - 1] = nodes[k];
				values[b - 1] = valuesAll[k]; 
			}
			else{
				sigma[b] = nodes[k];
				values[b] = valuesAll[k]; 
                //printVector(sigma);

			}
            //printVector(sigma);

			return 1;
		}
    }
}


double CalcPolynom(const std::vector<double>& coeffs, double x, int N){
    double ans = coeffs[1];
    for (int i = 1; i < N-1; ++i){
        double xi = pow(x, i);
        //std::cout<<coeffs[i+1];
        ans += xi * coeffs[i+1];
    }

    return ans; 
    
}



void  WriteToFile(double a, double b, const std::string& filename, std::vector<double>& coeffs){
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        
        int N = coeffs.size()-1;
        double h = (b - a) / (double)(3*(N));
        //h /= 3.;
        //int M = (double)(N-1)/h;
        double xi = a;
        double tmp1,tmp2,err;

        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x)"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"Pn"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"err"<<std::setw(9)<<" "<<std::endl;
        
        
        for (int i = 0; i < N-1; ++i){ 
            xi += h;  
            tmp1 = f(xi);
            tmp2 =  CalcPolynom(coeffs,xi,N);
            err = (tmp1-tmp2);
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;
            xi += h;   
            tmp1 = f(xi);
            tmp2 = CalcPolynom(coeffs,xi,N);
            err = (tmp1-tmp2);
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl; 
            xi += h;   
            tmp1 = f(xi);
            tmp2 = CalcPolynom(coeffs,xi,N);
            err = (tmp1-tmp2);
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;  

        }
        xi += h;  
            tmp1 = f(xi);
            tmp2 =  CalcPolynom(coeffs,xi,N);
            err = (tmp1-tmp2);
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

    }
    else {
        std::cerr << "Error opening file" << std::endl;
    }
}
