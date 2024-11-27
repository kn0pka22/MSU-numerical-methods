
#include "task_6.hpp"

double fun0(double x, double y){ return (x + y);}
double fun1(double x, double y){ return x * x + y * y;}
double fun2(double x, double y){ return (x * x * x * x + x * x * y * y + y * y * y *y);}
double fun3(double x, double y){ return x - y;}
double fun4(double x, double y){ return x * y;}
  


double FunctionWithName::TrueRes(double xa, double xb, double ya, double yb){
    double trueRes = 0;
    switch(id) {
        case 0: // f(x, y) = x + y
            trueRes = (xb*xb*yb + xb*yb*yb - xa*xa*ya - xa*ya*ya) / 2.;
            break;
        
    //     case 2: // f(x, y) = x^2 + y^2
    //         trueRes = (x2*x2*x2*y2 + x2*y2*y2*y2 - x1*x1*x1*y1 - x1*y1*y1*y1 ) / 3.;
    //         break;

    //     case 3: // f(x, y) = x^4 + x^2*y^2 + y^4
    //         trueRes =  x2*y2 * (pow(x2, 4) + pow(y2, 4)) / 5. + pow(x2, 3)* pow(y2, 3) / 9. - x1*y1 * (pow(x1, 4) + pow(y1, 4)) / 5. - pow(x1, 3)* pow(y1, 3) / 9. ;
    //         break;

    //     case 4: // f(x, y) = x - y
    //         trueRes =  (x2*x2*y2 - x2*y2*y2 - x1*x1*y1 + x1*y1*y1 ) / 2.;
    //         break;

        default: // f(x, y) = x * y
            trueRes = (xb * xb * yb * yb - xa * xa * ya * ya) / 4.;
            break;
    }

    return trueRes;
}

// void PointFromNum(int k, double& x, double& y, int N, double ){
//     x
// }

void triangulation(int N, double Lx, double Ly, const std::string &filename) {
    double hx = Lx / N;
    double hy = Ly / N;
    
    int number = 1;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(15);
    
    //out << "grid step sizes:       " << hx << " " << hy << "\n";

    out << "num of nodes:          " << (N + 1) * (N + 1) << "\n";
    // Printing the node coordinates
    for (int i = 0; i < (N + 1) * (N + 1); ++i) {
        out <<std::setw(4)<< i + 1 << ": " 
            <<std::setw(20)<< (i % (N + 1)) * hx << " " 
            <<std::setw(20)<< (i / (N + 1)) * hy << "\n";
    }
    
    
    out << "num of triangles:      " << 2 * N * N  << "\n";
    int LeftUpPoint = 0;
    int LeftDownPoint = 0;
    int RightDownPoint = 0;
    int RightUpPoint = 0;

    // Printing the triangles
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            

            LeftUpPoint    = i    *(N+1) + (j) + 1;
            RightUpPoint   = i    *(N+1) + (j + 1) + 1;
            LeftDownPoint  = (i+1)*(N+1) + j + 1;
            RightDownPoint = (i+1)*(N+1) + (j + 1) + 1;

            out << LeftUpPoint   << " "
                << RightUpPoint  << " "
                << LeftDownPoint << "\n";
            out << LeftDownPoint   << " "
                << RightUpPoint  << " "
                << RightDownPoint << "\n";

        }
    } 
    
    // out << "num of triangle sides: " << 3 * N * N + 2 * N << "\n";
    // // Printing the triangle edges
    // for (int j = 1; j < N + 1; ++j) {
    //     for (int i = 1; i < N + 1; ++i) {
    //         LeftUpPoint = i + (j - 1) * (N + 1);
            
    //         if (i == 1 && j == 1) {
    //             out << number << ": " << 2 << " " << N + 2 << "\n";
    //             number++;
    //         }
    //         else if (i == 1) {
    //             out << number << ": " << LeftUpPoint << " " << LeftUpPoint + 1 << "\n";
    //             number++;
    //             out << number << ": " << LeftUpPoint + 1 << " " << LeftUpPoint + N + 1 << "\n";
    //             number++;
    //         }
    //         else if (j == 1) {
    //             out << number << ": " << LeftUpPoint << " " << LeftUpPoint + N + 1 << "\n";
    //             number++;
    //             out << number << ": " << LeftUpPoint + 1 << " " << LeftUpPoint + N + 1 << "\n";
    //             number++;
    //         }
    //         else {
    //             out << number << ": " << LeftUpPoint << " " << LeftUpPoint + 1 << "\n";
    //             number++;
    //             out << number << ": " << LeftUpPoint + 1 << " " << LeftUpPoint + N + 1 << "\n";
    //             number++;
    //             out << number << ": " << LeftUpPoint + 1 << " " << LeftUpPoint + N + 1 << "\n";
    //             number++;
    //         }
    //     }
    // }

    
    // out << "num of edge nodes:     " << 4 * N << "\n\n";
    
   
    

   
    
    // number = 1;
    
    // // Printing the edge nodes
    // for (int j = 1; j < N + 1; ++j) {
    //     out << number << ": " << j << " " << j + 1 << "\n";
    //     number++;
    // }
    
    // for (int j = 0; j < N; ++j) {
    //     out << number << ": " << j * (N + 1) + 1 << " " << j * (N + 1) + N + 2 << "\n";
    //     number++;
    //     out << number << ": " << j * (N + 1) + 1 + N << " " 
    //         << j * (N + 1) + 1 + N + N + 1 << "\n";
    //     number++;
    // }
    
    // for (int j = 1; j < N + 1; ++j) {
    //     out << number << ": " << (N + 1) * N + j << " " 
    //         << (N + 1) * N + 1 + j << "\n";
    //     number++;
    // }
    
    out.close();
}



double IntegrateQuadr2(int N, std::function<double(double, double)> f) {
    double ans = 0.0;
    double anx, any, bnx, bny, cnx, cny;
    double h = 1.0 / double(N);

    for (int i = 1; i < N+1; ++i){  
        for (int j = 1; j < N+1; ++j){  
            anx = (i - 1) * h;
            any = (j - 1) * h;
            bnx = anx + h;
            bny = any;
            cnx = anx;
            cny = any + h;

            
            ans += 0.5 * h * h * (f((anx + bnx) / 2.0, (any + bny) / 2.0) +
                                  f((anx + cnx) / 2.0, (any + cny) / 2.0) +
                                  f((cnx + bnx) / 2.0, (cny + bny) / 2.0)) / 3.0;

            // The second triangle in the cell
            anx = (i) * h;
            any = (j) * h;
            bnx = anx - h;
            bny = any;
            cnx = anx;
            cny = any - h;

            ans += 0.5 * h * h * (f((anx + bnx) / 2.0, (any + bny) / 2.0) +
                                  f((anx + cnx) / 2.0, (any + cny) / 2.0) +
                                  f((cnx + bnx) / 2.0, (cny + bny) / 2.0)) / 3.0;
        }
    }

    return ans;
}





double IntegrateQuadr1(int N, std::function<double(double, double)> f) {
    double ans = 0.0;
    double anx, any, bnx, bny, cnx, cny;
    double h = 1.0 / double(N);
    double x, y;

    for (int i = 1; i < N+1; ++i){  
        for (int j = 1; j < N+1; ++j){  
            anx = (i - 1) * h;
            any = (j - 1) * h;
            bnx = anx + h;
            bny = any;
            cnx = anx;
            cny = any + h;

            x = (anx + bnx + cnx) / 3.0;
            y = (any + bny + cny) / 3.0;

            ans += 0.5 * h * h * f(x,y);

            // The second triangle in the cell
            anx = (i) * h;
            any = (j) * h;
            bnx = anx - h;
            bny = any;
            cnx = anx;
            cny = any - h;

            x = (anx + bnx + cnx) / 3.0;
            y = (any + bny + cny) / 3.0;

            ans += 0.5 * h * h * f(x,y);

        }
    }

    return ans;
}

