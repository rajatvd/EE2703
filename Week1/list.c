#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int main(){
    double n[1000];
    n[0]=0.2;
    double alpha = M_PI, temp;
    int k;
    for(k=1;k<1000;k++){
        temp = (n[k-1]+M_PI)*100;
        n[k] = temp - (long)(temp);
    }
    for(k=0;k<1000;k++){
        printf("%0.4f\n",n[k]);
    }
}