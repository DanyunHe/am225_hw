#include <cmath>
#include <cstdio>
#include <cstring>
#include "omp.h"

/** Custom random number generator based of "Ran" routine in Numerical Recipes
 * by Press et al. */
class custom_rng{
    public:
        unsigned long a,b,c; 
        custom_rng(unsigned long seed): b(4101842887655102017L), c(1) {
            if(sizeof(unsigned long)<8) {
                fputs("Error: 'unsigned long' type too short\n",stderr);
                exit(1);
            }
            a=seed^b;int64();
            b=a;int64();
            c=b;int64();
        }

        unsigned long int64() {
            a=a*2862933555777941757L+7046029254386353087L;
            b^=b>>17;b^=b<<31;b^=b>>8;
            c=4294957665U*(c&0xffffffff)+(c>>32);
            unsigned long d=a^(a<<21);
            d^=d>>35;d^=d<<4;
            return (d+b)^c;
        }
        inline double doub(){
            return 5.42101086242752217E-20*int64();
        }
};

const int m=80;
const int n=40;
// const int me=m+2,ne=n+2;
int* c=new int[m*n]; //cellular automaton on m by n grid 
int* cc=new int[m*n]; // next step cell grid 
custom_rng rng=custom_rng(0); // random number generator

// Initialize cell grid at time 0
void init(){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			if((i>=m/2-6)&&(i<m/2+6)&&(j>=n/2-6)&&(j<n/2+6)){
				if(rng.doub()<0.75){c[j+n*i]=1;}
				else{c[j+n*i]=0;}
			}
			else{c[j+n*i]=0;}
		}
	}
}

// Perform a gerneation step
void generate(){

	int nij;

#pragma omp parallel for 
	for(int i=0;i<m;i++){
		// Boundary condition
		// int il=(i==0)?m-1:i-1;
		// int ir=(i==m-1)?0:i+1;
		for(int j=0;j<n;j++){
			int il=(i==0)?m-1:i-1;
			int ir=(i==m-1)?0:i+1;
			int jl=(j==0)?n-1:j-1;
			int jr=(j==n-1)?0:j+1;

			// Calculate Nij for cell at (i,j)
			nij=c[ir*n+jl]+c[ir*n+j]+c[ir*n+jr]+c[i*n+jl]
				+c[i*n+jr]+c[il*n+jl]+c[il*n+j]+c[il*n+jr];
			int idx=j+i*n;
			// nij=c[idx-ne-1]+c[idx-ne]+c[idx-ne+1]+c[idx-1]
			// 	+c[idx+1]+c[idx+ne-1]+c[idx+ne]+c[idx+ne+1];

			// Update cell status 
			if(c[idx]==1){
				if(nij>=1&&nij<=5){cc[idx]=1;}
				else{cc[idx]=0;}
			}
			else if(nij==3){cc[idx]=1;}
			else {cc[idx]=0;}
		}
	}

	memcpy(c,cc,m*n*sizeof(int));
}

// Save the cell information 
void output(char* filename){

	FILE *fp=fopen(filename,"wb");
	if(fp==NULL){
		fputs("error: fail to open file",stderr);
		exit(-1);
	}

	fwrite(c,sizeof(int),n*m,fp);
	fclose(fp);

}

int main(){

	double t0=omp_get_wtime();
	init();
	int duration=151;
	for(int t=0;t<duration;t++){
		if(t%25==0){
			char fn[10];
			sprintf(fn,"%d.out",t);
			output(fn);
		}

		generate();
	}

	delete [] c;
	delete [] cc;

	double t1=omp_get_wtime();
	int num_thread=omp_get_max_threads();
	printf("running time for n %d,number of threads %d is %f\n",n,num_thread,(t1-t0)/duration);

	return 0;

}