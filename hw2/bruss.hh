#ifndef BRUSS_HH
#define BRUSS_HH

#include "sol_rk4d.hh"
#include "fsal.hh"

/** This class has functions to specify the test Brusselator problem. */
class bruss {
    public:
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void brus_ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void brus_init(double *q) {
            *q=1.5;
            q[1]=3.;
        }
};

class brus_rk4d : public rk4d, public bruss{
    public:
        brus_rk4d() : rk4d(2) {}
        ~brus_rk4d() {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};


class brus_fsal : public fsal, public bruss {
    public:
        brus_fsal() : fsal(2) {}
        ~brus_fsal() {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

#endif
