#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fsal.hh"
#include "sol_rk4d.hh"
#include "bruss.hh"
#include "omp.h"

int main() {
    // // Compute reference solution
    brus_rk4d br;
    br.solve_fixed(20.,20000000,false,0);
    double ref0=br.q[0],ref1=br.q[1];
    
    int fcount[101];
    double err[101];
    for(int i=0;i<=100;i++) {
        double lambda= 1e-3*pow(10.,-0.1*i);
        // int steps=int(100.*pow(1000.,0.01*i));
        brus_fsal sol_fsal; 
        sol_fsal.solve(20.,lambda,2000,0,false);
        double dy0=ref0-sol_fsal.q[0],
               dy1=ref1-sol_fsal.q[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=sol_fsal.fcount;
    }
    for(int i=0;i<=100;i++) printf("%d %g\n",fcount[i],err[i]);
    puts("\n");
}