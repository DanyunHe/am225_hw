#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <cmath>

#include "rbf.hh"
#include "omp.h"

int main() {

    for(int k=10;k<10000;k+=k>>2) {
        rbf r(k,2);
        r.set_length_scale(5./sqrt(k));
        r.init_random(0);
        int tk;
        double pk=0.;

        tk=r.assemble_matrix();

        double t0=omp_get_wtime(),t1,t_lapack,t_cg=0;
        int l=0;
        // Count the total non-zero entries
        
        do {
            r.solve_weights_lapack();
            l++;t1=omp_get_wtime();
        } while(t1<t0+0.5);
        t_lapack=(t1-t0)/l;

        r.make_table();
        t0=omp_get_wtime();l=0;
        int bls=sqrt(k);
        do {
            pk+=r.solve_weights_conj_grad(bls);
            l++;t1=omp_get_wtime();
        } while(t1<t0+1);
        t_cg=(t1-t0)/l;
        pk/=l;

        printf("%d %g %g %d %g\n",k,t_lapack,t_cg,tk,pk);
    }
}