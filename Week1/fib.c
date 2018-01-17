#include<stdlib.h>
#include<stdio.h>

int main(){
    int n=1, nold=1, new=0,k;
    printf("1 %d\n",n);
    printf("2 %d\n",nold);
    for(k=3;k<=10;k++){
        new = n+nold;
        nold=n;
        n=new;
        printf("%d %d\n",k,new);
    }
}