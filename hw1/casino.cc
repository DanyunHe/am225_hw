#include <cmath>
#include <cstdio>
#include <sys/types.h>
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

int main(){
	// Play 1e9 games and calculate the total winnings 
	unsigned int w=0;
	int n=1e9;
	double t0=omp_get_wtime();

	int num_thread=omp_get_max_threads();
	printf("number of threads %d\n",num_thread);
    // Create an array of pointers to random number generators
    custom_rng* c[num_thread];
    // Create random number generators, each with a different initial seed
    for(int i=0;i<num_thread;i++) c[i]=new custom_rng(i);

#pragma omp parallel for reduction(+:w)
	for(int i=0;i<n;i++){
		int tid=omp_get_thread_num();
		double x=0.;
		int win=0;
		while(x<1.){
			x+=c[tid]->doub();
			win++;
		}
		w+=win;
	}
	printf("average win is %f\n",w/1e7);
	
	for(int i=num_thread-1;i>=0;i--) delete [] c[i];

	double t1=omp_get_wtime();
	printf("running time %f\n",t1-t0);

}
