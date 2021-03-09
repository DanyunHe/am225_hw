#include "fasl.hh"
#include <cstdio>
#include <cstring>

/** Initialized the fsal solver.
 * \param[in] dof_ the number of degree of freedoms. */
fsal::fsal(int dof_):dof(dof_),fcount(0),t(0.),q(new double[dof]),
	dq(new double[dof]),
	k1(new double[dof]),k2(new double[dof]),k3(new double[dof]),k4(new double[dof]){}

/** The class destructor. */
fsal::~fsal(){
	delete [] k4;
	delete [] k3;
	delete [] k2;
	delete [] k1;
	delete [] dq;
}



/** Print the solution at given state. */
void fsal::print(double t_,double *in){
	printf("%g",t_);
	for(int i=0;i<dof;i++) printf("%g",in[i]);
}

void fsal::solve(){
	// Set up initial condition and compute timestep. Use the t_den variable to
    // mark where the dense output has got to.
    init();
    double dt=0.01,t_den=0,dt_den;

    // Perform integration steps
    if(output) print(t,q);
    if(d_steps>0) {
        dt_den=duration/d_steps;
        print(t,q);
    }

    ff(t,q,k1);
    for(int i=0;i<steps;i++) {
        
        dt=step(dt);

        // Do any dense output interpolation
        if(d_steps>0) {
            while(t_den+dt_den<t) {
                t_den+=dt_den;
                dense_output(1.+(t_den-t)/dt,dt);
                print(t_den,k3);
            }
        }

        // Move the required data for the next step into position
        memcpy(q,dq,dof*sizeof(double));
        memcpy(k1,k5,dof*sizeof(double));

        if(output) print(t,q);
    }

}

/** Computes a Hermite interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] theta the fraction of the timestep at which to evaluate the
 *                  interpolation.
 * \param[in] dt the length of the current timestep. */
void fsal::dense_output(double theta,double dt) {
    double mth=1-theta;

    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k1, and the new derivative is in
    // k2
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+theta*dq[i]
             -theta*mth*((1-2*theta)*(dq[i]-q[i])+dt*(theta*k2[i]-mth*k1[i]));
}

/** Step forward for a given step size dt. */
int fsal::step(double dt){

	// Second step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k1[i]/3.;
    ff(t+dt/3.,dq,k2);

    // Third step
    for(int i=0;i<dof;i++) dq[i]=q[i]-dt*k1[i]/3.+dt*k2[i];
    ff(t+2.*dt/3.,dq,k3);

    // Fourth step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]-k2[i]+k3[i]);
    ff(t+dt,dq,k4);

	// Fifth step 
	fcount+=4;
	for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i])/8.;
    ff(t+dt,dq,k5);

    // Find adaptive step 
    while(true){

    }
    

    // Complete solution
    t+=dt;
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(1/12.)*(k1[i]+6*k2[i]+3*k3[i]+2*k5[i]);

    // Reuse k2 to store the derivative at the new solution
    // ff(t,dq,k2);

}

/** Find the adaptive step size. */
double fsal::step_size(){



	

}





