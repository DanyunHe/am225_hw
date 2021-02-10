#include <cmath>
#include <cstdio>

/** Find all numbers less than n
 * \param[in] A A[i]=1 implies (i+2) is prime. */
void sieve(int n,int* A){
	for(int i=0;i<n-1;i++){A[i]=1;}
	for(int i=0;i<int(sqrt(n))-1;i++){
		if(A[i]){
			int i2=(i+2)*(i+2);
			for(int j=i2;j<n;j+=(i+2)){
				A[j-2]=0;
			}
		}
	}
}

/** Perform division with remainder
 * \param[in] d,K The dk sequence for M.
 * \param[in] B the base.
 * \param[in] p the divisor.
 * \param[in] e the return solution for M div p.
 * \param[in] a the remainder M%p. */
void division(int* d,int K,int B,int p,int* e,int& a){
	e[K]=d[K]/p;
	a=d[K]%p;
	for(int k=K-1,k>=0;k--){
		alpha=a*B+d[k];
		e[k]=alpha/p;
		a=alpha%p;
	}	
}

int count(int* d,int K,int B,int n=2e5){

	int* A=new int[n-1];
	// Find all prime numbers less than n
	sieve(n,A);
	int count=0;
	for(int i=0;i<n-1;i++){
		if(A[i]==1){
			int* e=new int[K+1];
			int a;
			division(d,K,B,i+2,e,&a);
			if(a==0) count++;
			delete [] e;
		}
	}
	delete [] A;
}




int main(){

	int n=2*1e5;
	int* A=new int[n-1];
	// Find all prime numbers less than n
	sieve(n,A);


}