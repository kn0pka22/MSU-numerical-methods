#include <iostream>


int IntSum2n(){
    int S=0;
    int Sn=0;
    int qn=1;
    while(1){
        Sn = S + qn;
        if (S == Sn) break;
        qn = qn*2;
        S = Sn;
    }
    return S;
}

int deg(int x, int n){
    int res = 1; 
    for (int i=0;i<n;++i){
        res *= x; 
    }
    return res;
}

int func(int x){
    int degg = 1;
    //(x-2)^9
    
    return (deg(x, 9) - 9*deg(x,8)*deg(2,1) + 36*deg(x, 7) + 84*deg(x,6)+126*deg)
} 



int main(){
    std::cout<<IntSum2n()<<std::endl;

    return 0;
}
