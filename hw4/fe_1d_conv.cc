#include "cubic_1d_fe.hh"
#include "alter_1d_fe.hh"
#include "omp.h"

#include <cmath>

int main() {

    // Loop over a variety of different grid sizes
#pragma omp parallel for schedule(dynamic) ordered
    for(int i=0;i<=30;i++) {

        // Create the finite-element problem and set the Neumann boundary
        // condition to match the manufactured solution
        int j=int(10*pow(100,(1/30.)*i)+0.5);
        
        double t0=omp_get_wtime(),t1;
        alter_1d_fe cf(j);
        cf.g=exp(-1)*5*M_PI;
        // Initialize the source term for the manufactured solution, solve, and
        // the print the L2 error
        cf.init_mms();
        int count=0;
        do{
            cf.solve();
            t1=omp_get_wtime();
            count++;
        }while(t1-t0<0.5);
        

#pragma omp ordered
        printf("%d %g %g %g\n",j,cf.h,cf.l2_norm_mms(),(t1-t0)/count);
    }
}
