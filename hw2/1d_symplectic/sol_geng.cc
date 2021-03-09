#include "sol_geng.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Initializes the fourth-order Hammer--Hollingsworth solver, allocating
 * memory and setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
geng::geng(int dof_) : sol_base(dof_), dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k1b(new double[dof]), k2b(new double[dof]), k3b(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
geng::~geng() {
    delete [] k3b;
    delete [] k2b;
    delete [] k1b;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
}

/** Performs an integration step with the fourth-order Hammer--Hollingsworth solver.
 * \param[in] dt the integration step. */
bool geng::step(double dt) {
    int iter=0;
    double delsq,d,*c;
    const double r6=sqrt(6.);

    // Clear steps
    for(int i=0;i<dof;i++) k1[i]=k2[i]=k3[i]=0;

    do {

        // Check for too many iterations
        if(++iter>1000) {
            fputs("Too many iterations in IRK\n",stderr);
            exit(1);
        }

        // Perform update
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*((16.-r6)/72.*k1[i]+(328.-167.*r6)/1800.*k2[i]+(-2.+3.*r6)/450.*k3[i]);
        ff(t+dt*(4.-r6)/10.,dq,k1b);
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*((328.+167.*r6)/1800.*k1[i]+(16.+r6)/72.*k2[i]+(-2.-3.*r6)/450.*k3[i]);
        ff(t+dt*(4.+r6)/10.,dq,k2b);
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*((85.-10.*r6)/180.*k1[i]+(85.+10.*r6)/180.*k2[i]+1./18.*k3[i]);
        ff(t+dt,dq,k3b);
        fcount+=3;

        // Find size of step from previous iteration
        delsq=0;
        for(int i=0;i<dof;i++) {
            d=k1[i]-k1b[i];delsq+=d*d;
            d=k2[i]-k2b[i];delsq+=d*d;
            d=k3[i]-k3b[i];delsq+=d*d;
        }

        // Switch k1<->k1b and k2<->k2b array pointers. This will make k1 & k2
        // be used as the new values on the next iteration.
        c=k1b;k1b=k1;k1=c;
        c=k2b;k2b=k2;k2=c;
        c=k3b;k3b=k3;k3=c;
    } while(delsq>1e-25);

    // Complete solution
    for(int i=0;i<dof;i++) q[i]+=dt*((16.-r6)/36.*k1[i]+(16.+r6)/36.*k2[i]+1./9.*k3[i]);
    t+=dt;
    return true;
}
