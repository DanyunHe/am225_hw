#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "fsal.hh"
#include "sol_rk4d.hh"
#include "disc.hh"


int main() {
    
    double lambda=3./1000;
    double duration=48.+exp(-1.);

    int fcount[101];
    double err[101];

        // int steps=int(100.*pow(1000.,0.01*i));
    // disc_rk4d dr;
    int step;
    for(int i=0;i<101;i++){
        disc_rk4d dr;
        step=int(1000.*pow(10000.,0.01*i));
        dr.solve_fixed(duration,step,false,0);
        double dy0=1.-dr.q[0],
               dy1=exp(-1.)-dr.q[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=dr.fcount;

    }
    for(int i=0;i<101;i++) printf("%d %g\n",fcount[i],err[i]);
    puts("\n");
    // printf("%d %g\n",fcount,err);
    // puts("\n");
}


/**
int main() {
    // // Compute reference solution
    // brus_rk4d br;
    // br.solve_fixed(20.,20000000,false,0);
    // double ref0=br.q[0],ref1=br.q[1];
 
    // double lambda=3./1000;
    double duration=48.+1./2.718281828459045;
    double ref0=1.,ref1=1./2.718281828459045;

    int fcount[101];
    double err[101];

        // int steps=int(100.*pow(1000.,0.01*i));
    // disc_fsal dr;
    int step;
    double lambda;
    for(int i=0;i<101;i++){
        disc_fsal dr;
        // printf("i %d, 0.1i %g\n",i,-0.1*i);
        lambda=1e-2*pow(10.,-0.1*i);
        dr.solve(duration,lambda,100000,0,false);
        double dy0=ref0-dr.q[0],
               dy1=ref1-dr.q[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=dr.fcount;

    }
    for(int i=0;i<101;i++) printf("%d %g\n",fcount[i],err[i]);
    puts("\n");
    // printf("%d %g\n",fcount,err);
    // puts("\n");
}
*/


