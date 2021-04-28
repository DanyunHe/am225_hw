#include "bilinear_2d_fe.hh"

#include <cmath>

int main() {

    // Loop over a variety of different grid sizes
#pragma omp parallel for schedule(dynamic) ordered
    for(int i=0;i<=20;i++) {

        // Create the finite-element problem 
        int j=int(10*pow(100,(1/60.)*i)+0.5);

        // Construct the finite-element class
        bilinear_2d_fe bf(j);
        
        // Initialize the source function 
        bf.init();

        // Solve the finite-element problem using the conjugate gradient method
        bf.solve();
        // bf.print();
#pragma omp ordered
        printf("%d %g %g\n",j,bf.h,bf.l2_norm_mms());
    }
}
