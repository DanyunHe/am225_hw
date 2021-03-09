#include "orbit.hh"

#include <cstdio>
#include <cmath>

/** Sets up the initial conditions for the ODE.
 * \param[in] q the array to write to. */
void orbit::orb_init(double *q) {
    // const double a=1.;
    // const double e=1/3;
    // const double r=a*(1+e);
    // q 0 1 2 are initial velocity
    // q 3 4 5 are initial position
    *q=0; // p1
    q[1]=1.68883701;  // p2
    q[2]=0.2;  // p3
    q[3]=2.5; // q1
    q[4]=0; // q2
    q[5]=0; // q3
    init_h=hamiltonian(q);
    //printf("# Initial H(p,q)=%.12g\n",init_h);
}

/** Prints the current state of the solution, plus the Hamiltonian.
 * \param[in] t_ the current simulation time.
 * \param[in] q the solution array. */
void orbit::orb_print(double t_,double *q) {
    printf("%g %g %g %g %g %g %g %.12g\n",t_,*q,q[1],q[2],q[3],q[4],q[5],hamiltonian(q)-init_h);
}

/** Computes the Hamiltonian.
 * \param[in] q the solution array.
 * \return The Hamiltonian. */
double orbit::hamiltonian(double *q) {
    const double Omega = 0.25;
    const double A = 1.0;
    const double C = 1.0;
    const double a = 1.25;
    const double b = 1.0;
    const double c = 0.75;
    return 0.5*(*q*(*q)+q[1]*q[1]+q[2]*q[2])+Omega*(*q*q[4]-q[1]*q[3])+A*log(C+q[3]*q[3]/(a*a)+q[4]*q[4]/(b*b)+q[5]*q[5]/(c*c));
}
