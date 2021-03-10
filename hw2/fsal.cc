#include "fsal.hh"
#include <cstdio>
#include <cstring>
#include <cmath>

/** Initialized the fsal solver.
 * \param[in] dof_ the number of degree of freedoms. */
fsal::fsal(int dof_):dof(dof_),fcount(0),t(0.),
    fac(0.9),facmax(3.),facmin(1./3.),atol(0.),rtol(0.),
    q(new double[dof]),dq(new double[dof]),
	k1(new double[dof]),k2(new double[dof]),
    k3(new double[dof]),k4(new double[dof]),k5(new double[dof]){}

/** The class destructor. */
fsal::~fsal(){
    delete [] k5;
	delete [] k4;
	delete [] k3;
	delete [] k2;
	delete [] k1;
	delete [] dq;
    delete [] q;
}



/** Print the solution at given state. */
void fsal::print(double t_,double *in){
	printf("%g",t_);
	for(int i=0;i<dof;i++) printf(" %g",in[i]);
    puts("");
}

double fsal::absolute(double x){
    if(x>=0) return x;
    else return -x;
}
/** Solve ODE system with adaptive step.
 * \param[in] T the duration.
 * \param[in] lambda the value for Atol and Rtol.
 * \param[in] n number of steps.
 * \param[in] dn number of dense output steps.
 * \param[in] output whether print each integration step. */
void fsal::solve(double T,double lambda,int n,int dn,bool output){
	// Set up initial condition and compute timestep. Use the t_den variable to
    // mark where the dense output has got to.
    init();
    double dt=T/n,t_den=0,dt_den;

    // Perform integration steps
    if(output) print(t,q);
    if(dn>0) {
        dt_den=T/dn;
        print(t,q);
    }

    ff(t,q,k1);
    bool last=false;
    while(t<T){
        // Check if it is the last step
        double old_T = t;
        if(t<T-dt){
            dt=step(dt,lambda,last);
        }
        else{
            dt=T-t;
            last=true;
            dt=step(dt,lambda,last);
        }

        double dt_taken = t - old_T;

        // Do any dense output interpolation
        if(dn>0) {
            while(t_den+dt_den<t) {
                t_den+=dt_den;
                dense_output(1.+(t_den-t)/dt_taken,dt_taken);
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
    // k5
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+theta*dq[i]
             -theta*mth*((1-2*theta)*(dq[i]-q[i])+dt*(theta*k5[i]-mth*k1[i]));
}

/** Step forward for a given step size dt. */
double fsal::step(double dt,double lambda,bool last){

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
    

    // Compute error
    double err=0.,qq=3.;
    double yi,yhat,sc;
    atol=lambda;
    rtol=0;
    for(int i=0;i<dof;i++){
        yi=dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i])/8.;
        yhat=dt*(1/12.)*(k1[i]+6*k2[i]+3*k3[i]+2*k5[i]);
        sc=atol+rtol*fmax(absolute(q[i]),absolute(q[i]+yi));
        err+=(yi-yhat)*(yi-yhat)/(sc*sc);
    }
    err=pow((err/dof),0.5);
    

    // Compute an adaptive step size
    while(err>1&&!last){
        // reject
        dt=dt*fmin(facmax,fmax(facmin,fac*pow((1./err),1./(qq+1.))));
        err=0;

        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k1[i]/3.;
        ff(t+dt/3.,dq,k2);

        // Third step
        for(int i=0;i<dof;i++) dq[i]=q[i]-dt*k1[i]/3.+dt*k2[i];
        ff(t+2.*dt/3.,dq,k3);

        // Fourth step
        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]-k2[i]+k3[i]);
        ff(t+dt,dq,k4);

        for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i])/8.;
        ff(t+dt,dq,k5);
        fcount+=4;

        double yi,yhat,sc;
        for(int i=0;i<dof;i++){
            yi=dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i])/8.;
            yhat=dt*(1/12.)*(k1[i]+6*k2[i]+3*k3[i]+2*k5[i]);
            sc=atol+rtol*fmax(absolute(q[i]),absolute(q[i]+yi));
            err+=(yi-yhat)*(yi-yhat)/(sc*sc);
        }
        err=pow((err/dof),0.5);

    }
    

    // Complete solution
    t+=dt;
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i])/8.;
    dt=dt*fmin(facmax,fmax(facmin,fac*pow((1./err),1./(qq+1.))));

    // printf("err %g dt %g lambda %g\n",err,dt,lambda);

    return dt;

}








