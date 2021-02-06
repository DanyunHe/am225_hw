#include <cmath>
#include <cstdio>

void sieve(int n){
	int A[n-1];
	for(int i=0;i<n-1;i++){A[i]=1;}
	for(int i=0;i<int(sqrt(n))-1;i++){
		if(A[i]){
			int i2=(i+2)*(i+2);
			for(int j=i2;j<n;j+=(i+2)){
				A[j-2]=0;
			}
		}
	}
	//all i, that A[i]=false => i+2 is prime

}

void devision(){
	
}