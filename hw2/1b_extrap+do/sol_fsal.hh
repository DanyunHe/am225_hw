#ifndef SOL_FSAL_HH
#define SOL_FSAL_HH

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class fsal {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        fsal(int dof_);
        virtual ~fsal();
        void print(double t_,double *in);
        void dense_output(double theta,double dt);
        void solve_fixed(double t_end,int iters,bool output,int d_steps, double lambda);
        double step(double dt, double lambda, int islast);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        double *dq;
        double *dq0;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        double *k5;
};

#endif
