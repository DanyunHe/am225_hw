#include <cmath>
#include <cstdio>
#include "omp.h"

/** Find all numbers less than n
 * \param[in] A A[i]=1 implies (i+2) is prime. */
void sieve(int n,int* A){
	for(int i=0;i<n-1;i++){A[i]=1;}
	for(int i=0;i<int(sqrt(n))-1;i++){
		if(A[i]){
			int i2=(i+2)*(i+2);
			for(int j=i2;j<=n;j+=(i+2)){
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
	int alpha;
	for(int k=K-1;k>=0;k--){
		alpha=a*B+d[k];
		e[k]=alpha/p;
		a=alpha%p;
	}	
}
/** Count number of primes less than n that are factors of M.
 * \param[in] d the sequence of dk for M.
 * \param[in] K the size of dk. */
int count(int* d,int K,int B,int n){

	int* A=new int[n-1];
	// Find all prime numbers less than n
	sieve(n,A);
	int result=0;
#pragma omp parallel for reduction(+:result) 
	for(int i=0;i<n-1;i++){
		if(A[i]==1){
			int* e=new int[K+1];
			int a;
			division(d,K,B,i+2,e,a);
			if(a==0) {
				result++;
				printf("the factor is %d\n",i+2);
			}
			delete [] e;
		}
	}
	delete [] A;
	return result;
}


int main(){

	// Initialize base B value
	int B_pow=10;
	int B=pow(2,B_pow);

	// N=2^n-1
	int n=82589933;

	// Convert N into the expression in Eq1
	// Calculate dk
	int a=n/B_pow+1,b=n%B_pow;
	int* d=new int[a];
	for(int i=0;i<a-1;i++){
		d[i]=B-1;
	}
	d[a-1]=pow(2,b)-1;

	// Count number of primes that are factors of N
	int result=count(d,a-1,B,2e5);
	printf("n %d B pow %d result %d\n",n,B_pow,result);


}