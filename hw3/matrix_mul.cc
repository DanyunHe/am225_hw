#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "blas.h"
#include "omp.h"


/** Standard matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void std_mul(int n,double* A, double* B, double* C){
	
	double temp;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			temp=0.;
			for(int k=0;k<n;k++){
				temp+=A[k*n+i]*B[j*n+k];
			}
			C[j*n+i]=temp;
		}
	}

}

/** Strassen's algorithm matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void strassen_mul(int n, double* A,double* B,double* C){

	if(n==1){
		C[0]=A[0]*B[0];
		return;
	}

	// Construct submatrix
	int nn=n/2,idx;
	const int nn2=nn*nn;
	double *A00=new double[nn2],*A01=new double[nn2],*A10=new double[nn2],*A11=new double[nn2];
	double *B00=new double[nn2],*B01=new double[nn2],*B10=new double[nn2],*B11=new double[nn2];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			idx=j*n+i;
			if(i<nn){
				if(j<nn){
					A00[j*nn+i]=A[idx];
					B00[j*nn+i]=B[idx];
				}
				else{
					A01[(j-nn)*nn+i]=A[idx];
					B01[(j-nn)*nn+i]=B[idx];
				}
			}
			else if(j<nn){
				A10[j*nn+i-nn]=A[idx];
				B10[j*nn+i-nn]=B[idx];

			}
			else{
				A11[(j-nn)*nn+i-nn]=A[idx];
				B11[(j-nn)*nn+i-nn]=B[idx];
			}
		}
	}

	// printf("%g %g %g %g\n",A00[0],A01[0],A10[0],A11[0]);

	// Calculate Q
	double *Q0=new double[nn2],*Q1=new double[nn2],*Q2=new double[nn2],*Q3=new double[nn2],\
		*Q4=new double[nn2],*Q5=new double[nn2],*Q6=new double[nn2];

	double *temp1=new double[nn2],*temp2=new double[nn2];
	
	// Q0
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A00[j*nn+i]+A11[j*nn+i];
			temp2[j*nn+i]=B00[j*nn+i]+B11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q0);

	// Q1
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A10[j*nn+i]+A11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,B00,Q1);

	// Q2
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=B01[j*nn+i]-B11[j*nn+i];
		}
	}
	strassen_mul(nn,A00,temp1,Q2);

	// Q3
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=B10[j*nn+i]-B00[j*nn+i];
		}
	}
	strassen_mul(nn,A11,temp1,Q3);


	// Q4
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A00[j*nn+i]+A01[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,B11,Q4);

	// Q5
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A10[j*nn+i]-A00[j*nn+i];
			temp2[j*nn+i]=B00[j*nn+i]+B01[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q5);

	// Q6
	for(int i=0;i<nn;i++){
		for(int j=0;j<nn;j++){
			temp1[j*nn+i]=A01[j*nn+i]-A11[j*nn+i];
			temp2[j*nn+i]=B10[j*nn+i]+B11[j*nn+i];
		}
	}
	strassen_mul(nn,temp1,temp2,Q6);
	// printf("%g %g %g %g %g %g %g\n",Q0[0],Q1[0],Q2[0],Q3[0],Q4[0],Q5[0],Q6[0]);


	// Calculate C
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i<nn){
				if(j<nn){
					idx=j*nn+i;
					C[j*n+i]=Q0[idx]+Q3[idx]-Q4[idx]+Q6[idx];
				}
				else{
					idx=(j-nn)*nn+i;
					C[j*n+i]=Q2[idx]+Q4[idx];
				}
			}
			else if(j<nn){
				idx=j*nn+i-nn;
				C[j*n+i]=Q1[idx]+Q3[idx];
			}
			else{
				idx=(j-nn)*nn+i-nn;
				C[j*n+i]=Q0[idx]+Q2[idx]-Q1[idx]+Q5[idx];
			}
		}
	}

	delete [] temp2;
	delete [] temp1;
	delete [] Q6;
	delete [] Q5;
	delete [] Q4;
	delete [] Q3;
	delete [] Q2;
	delete [] Q1;
	delete [] Q0;
	delete [] B11;
	delete [] B10;
	delete [] B01;
	delete [] B00;
	delete [] A11;
	delete [] A10;
	delete [] A01;
	delete [] A00;


}

/** BLAS matrix matrix multiplication.
 *\ param[in] n matrix dimension.
 *\ Calculate C=AB. */
void blas_mul(int n,double* A,double* B,double* C){
	char trans='n';
	double alpha=1.,beta=0.;

	dgemm_(&trans,&trans,&n,&n,&n,&alpha,A,&n,B,&n,&beta,C,&n);

}

/** Generate uniform random numbers in [-1,1]. */
double rnd(){
	return -1+2./RAND_MAX*static_cast<double>(rand());

}

/** Generate random matrix. */
void fill_matrix(int n,double* A){
	for(int i=0;i<n*n;i++){
		A[i]=10.*rnd();
	}

}

int main(){
	int n=8,l;
	double t0,t1,t_sr,t_std,t_blas;
	for(int i=0;i<7;i++){
		n*=2;
		// Generate random matrix
		double *A=new double[n*n],*B=new double[n*n],*C=new double[n*n];
		fill_matrix(n,A);
		fill_matrix(n,B);
		
		// Perform strassen's algorithm
		t0=omp_get_wtime();
		l=0;
		do{
			strassen_mul(n,A,B,C);
			l++;
			t1=omp_get_wtime();
		}while(t1<t0+0.5);
		t_sr=(t1-t0)/l;

		// Perform standard algorithm
		fill_matrix(n,A);
		fill_matrix(n,B);
		t0=omp_get_wtime();
		l=0;
		do{
			std_mul(n,A,B,C);
			l++;
			t1=omp_get_wtime();
		}while(t1<t0+0.5);
		t_std=(t1-t0)/l;

		// Perform BLAS
		fill_matrix(n,A);
		fill_matrix(n,B);
		char trans='n';
		double alpha=1.,beta=0.;
		t0=omp_get_wtime();
		l=0;
		do{
			fill_matrix(n,A);
			fill_matrix(n,B);
			// blas_mul(n,A,B,C);
			
			dgemm_(&trans,&trans,&n,&n,&n,&alpha,A,&n,B,&n,&beta,C,&n);
			l++;
			t1=omp_get_wtime();
		}while(t1<t0+0.5);
		t_blas=(t1-t0)/l;

		printf("%d %g %g %g\n",i,t_sr,t_std,t_blas);
		
		delete [] C;
		delete [] B;
		delete [] A;

	}
}







