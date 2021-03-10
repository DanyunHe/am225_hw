#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fsal.hh"
#include "sol_rk4d.hh"
#include "osc.hh"
#include "omp.h"

int main() {
    // // Compute reference solution
    // brus_rk4d br;
    // br.solve_fixed(20.,20000000,false,0);
    // double ref0=br.q[0],ref1=br.q[1];
    
    int fcount;
    double err;
    double lambda=3./1000;
    double duration=8.;

        // int steps=int(100.*pow(1000.,0.01*i));
    osc_fsal sol_fsal; 
    // sol_fsal.solve(duration,lambda,2000,1200,true);
    sol_fsal.solve(duration,lambda,2000,1200,true);
    double dy0=sol_fsal.sol0(duration)-sol_fsal.q[0],
            dy1=sol_fsal.sol1(duration)-sol_fsal.q[1];
    err=sqrt(dy0*dy0+dy1*dy1);
    fcount=sol_fsal.fcount;
    
    // printf("%d %g\n",fcount,err);
    // puts("\n");
}