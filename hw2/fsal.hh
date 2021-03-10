#ifndef FSAL_HH
#define FSAL_HH

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class fsal {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        double fac;
        double facmax;
        double facmin;
        double atol;
        double rtol;

        /** The solution vector. */
        double *q;
        fsal(int dof_);
        virtual ~fsal();
        double absolute(double x);
        void print(double t_,double *in);
        void dense_output(double theta,double dt);
        void solve(double T,double lambda,int n,int dn,bool output);
        double step(double dt,double lambda,bool last);
        virtual void init()=0;
        virtual void ff(double t_,double *in,double *out)=0;
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        double *k5;
};

#endif
