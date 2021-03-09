#ifndef SOL_GENG_HH
#define SOL_GENG_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using the fourth-order Hammer-Hollingsworth
 * method. */
class geng : public sol_base {
    public:
        geng(int dof_);
        virtual ~geng();
        virtual bool step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k1b;
        double *k2b;
        double *k3b;
};

#endif
