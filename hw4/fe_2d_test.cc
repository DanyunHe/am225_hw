#include <cmath>

#include "bilinear_2d_fe.hh"

int main() {

    // Construct the finite-element class
    bilinear_2d_fe bf(8);
    
    // Initialize the source function 
    bf.init();

    // Optional command to print the matrix in text form
    bf.print_matrix();

    // Solve the finite-element problem using the conjugate gradient method
    bf.solve();
    bf.print();
    printf("%g %g\n",bf.h,bf.l2_norm_mms());
}
