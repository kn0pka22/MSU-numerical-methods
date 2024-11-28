
#include "task_6.hpp"


Triangle::Triangle(int v1, int v2, int v3){
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
}

double PointFromNumVert(int numVert, int N, double h, char c){
    double res = (c == 'x') ? (((numVert-1) % (N + 1)) * h ) : (((numVert-1) / (N + 1)) * h);   
    //std::cout<<"c = "<<c<<" numVert = "<<numVert<<" N = "<<" res = "<<res<<std::endl;
    return res;
}


double FunctionWithName::TrueRes(double xa, double xb, double ya, double yb){
    double trueRes = 0;
    switch(id) {
        case 0: // f(x, y) = x + y
            trueRes = (yb-ya) * (xb-xa) * (xa + xb + ya + yb) / 2.;
            break;
        
        case 1: // f(x, y) = x^2 + y^2
            trueRes = (xb*xb*xb*yb + xb*yb*yb*yb - xa*xa*xa*ya - xa*ya*ya*ya ) / 3.;
            break;

        case 2: // f(x, y) = x^4 + x^2*y^2 + y^4
            trueRes =  xb*yb * (pow(xb, 4) + pow(yb, 4)) / 5. + pow(xb, 3)* pow(yb, 3) / 9. - xa*ya * (pow(xa, 4) + pow(ya, 4)) / 5. - pow(xa, 3)* pow(ya, 3) / 9. ;
            break;

        case 3: // f(x, y) = x - y
            trueRes =  (xb*xb*yb - xb*yb*yb - xa*xa*ya + xa*ya*ya ) / 2.;
            break;

        case 5:  // f(x, y) = sin(x * y)
            trueRes = -1. * sin(xb*yb)/(xb*yb) + sin(xa*ya)/(xa*ya);
;
            break;


        default: // f(x, y) = x - y
            trueRes = (xb * xb * yb * yb - xa * xa * ya * ya) / 4.;
            break;
    }

    return trueRes;
}

// void PointFromNum(int k, double& x, double& y, int N, double ){
//     x
// }



void triangulation(int N, double xa, double xb, double ya, double yb, const std::string &filename) {

    double hx = (xb-xa) / (double)N;
    double hy = (yb-ya) / (double)N;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }
    out << std::fixed << std::setprecision(15);
    
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

    Triangle tr1(LeftUpPoint, RightUpPoint, LeftDownPoint);
    Triangle tr2(RightUpPoint, RightDownPoint, LeftDownPoint);

    // Printing the triangles
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            

            tr1.v1 = i    *(N+1) + (j) + 1;
            tr1.v2 = i    *(N+1) + (j + 1) + 1;
            tr1.v3 = (i+1)*(N+1) + j + 1;

            tr2.v1 = i    *(N+1) + (j + 1) + 1;
            tr2.v2 = (i+1)*(N+1) + (j + 1) + 1;
            tr2.v3 = (i+1)*(N+1) + j + 1;

            out << tr1.v1 << " " << tr1.v2 << " " << tr1.v3 << "\n";
            out << tr2.v1 << " " << tr2.v2 << " " << tr2.v3 << "\n";
        }
    } 
    out.close();
}




double IntegrateQuadr1(int N, double xa, double xb, double ya, double yb, std::function<double(double, double)>& f, const std::string& filename){

    double res = 0.;

    double hx = (xb-xa) / (double)N;
    double hy = (yb-ya) / (double)N;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading!" << std::endl;
        return 1;
    }

    std::string line;
    bool readingTriangles = false;
    std::vector<Triangle> triangles;  

    while (getline(file, line)) {
        // We are looking for the line "num of triangles: .." or similar
        if (line.find("num of triangles:") != std::string::npos) {
            readingTriangles = true;  
            continue;  
        }

        if (readingTriangles) {
            std::stringstream ss(line);
            int v1, v2, v3;
            if (ss >> v1 >> v2 >> v3) {
                triangles.push_back(Triangle(v1, v2, v3));
            }
        }
    }

    file.close();  

    double x, y;
    for (int i = 0; i < triangles.size(); ++i){
        x = (PointFromNumVert(triangles[i].v1, N, hx, 'x')
          +  PointFromNumVert(triangles[i].v2, N, hx, 'x') 
          +  PointFromNumVert(triangles[i].v3, N, hx, 'x')) / 3.0;

        y = (PointFromNumVert(triangles[i].v1, N, hy, 'y')
          +  PointFromNumVert(triangles[i].v2, N, hy, 'y') 
          +  PointFromNumVert(triangles[i].v3, N, hy, 'y')) / 3.0;

        res += 0.5 * hx * hy * f(x,y);  
    }

    return res;
}




double GenereteFileForPCalculation(int numTests, double xa, double xb, double ya, double yb, 
                                FunctionWithName f, int N,
                                const std::string& fileForTriangulation,
                                const std::string& fileForP){



    std::ofstream fout(fileForP);
    if (!fout) {
        std::cerr << "Error opening file for p!" << std::endl;
        return -1;
    }


    fout << std::left 
         << std::setw(15) << "N" 
         << std::setw(20) << "Quadrature" 
         << std::setw(20) << "Analytical" 
         << std::setw(20) << "error" 
         << std::endl;

    double res1, res2, err1, err2;
    double tmp=0.;
    double p;

    std::vector<double> errVec;
    std::vector<double> NVec;

    for(int i = 0; i < numTests; ++i){

        triangulation(N, xa, xb, ya, yb, fileForTriangulation);
        res1 = IntegrateQuadr1(N, xa, xb, ya, yb, f.func, fileForTriangulation);
        res2 = f.TrueRes(xa, xb, ya, yb); 
        err1 = fabs(res1-res2);
        errVec.push_back(err1);
        NVec.push_back(N);

        fout << std::left 
             << std::setw(15) << std::fixed << std::setprecision(10) << N
             << std::setw(20) << std::fixed << std::setprecision(10) << res1
             << std::setw(20) << std::fixed << std::setprecision(10) << res2
             << std::setw(20) << std::fixed << std::setprecision(10) << err1
             << std::endl;

        N++; 

    }
    //std::cout<<N<<std::endl;

    fout.close();

    for (int i=1; i<errVec.size(); ++i){
        tmp += log(errVec[i-1]/errVec[i])/log(NVec[i]/NVec[i-1]); 
        //std::cout<<"tmp = "<< log(NVec[i]/NVec[i-1])<<std::endl;   
    }
    p = tmp / (double)(numTests-1);
    std::cout<<"p = "<<p<<std::endl;
    return p;
}


// double IntegrateQuadr2(int N, std::function<double(double, double)> f) {
//     double ans = 0.0;
//     double v1_x, v1_y, v2_x, v2_y, v3_x, v3_y;
//     double h = 1.0 / double(N);

//     for (int i = 1; i < N+1; ++i){  
//         for (int j = 1; j < N+1; ++j){  
//             v1_x = (i - 1) * h;
//             v1_y = (j - 1) * h;
//             v2_x = v1_x + h;
//             v2_y = v1_y;
//             v3_x = v1_x;
//             v3_y = v1_y + h;

            
//             ans += 0.5 * h * h * (f((v1_x + v2_x) / 2.0, (v1_y + v2_y) / 2.0) +
//                                   f((v1_x + v3_x) / 2.0, (v1_y + v3_y) / 2.0) +
//                                   f((v3_x + v2_x) / 2.0, (v3_y + v2_y) / 2.0)) / 3.0;

//             // The second triangle in the cell
//             v1_x = (i) * h;
//             v1_y = (j) * h;
//             v2_x = v1_x - h;
//             v2_y = v1_y;
//             v3_x = v1_x;
//             v3_y = v1_y - h;

//             ans += 0.5 * h * h * (f((v1_x + v2_x) / 2.0, (v1_y + v2_y) / 2.0) +
//                                   f((v1_x + v3_x) / 2.0, (v1_y + v3_y) / 2.0) +
//                                   f((v3_x + v2_x) / 2.0, (v3_y + v2_y) / 2.0)) / 3.0;
//         }
//     }

//     return ans;
// }


