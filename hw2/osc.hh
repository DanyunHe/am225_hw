#ifndef OSC_HH
#define OSC_HH

#include "sol_rk4d.hh"
#include "fsal.hh"

/** A 2-component oscillating ODE system. */
class osc {
    public:
        inline void osc_init(double *q) {
            *q=1;q[1]=0;
        }
        inline void osc_ff(double tt,double *in,double *out) {
            out[1]=tt*(*in);
            *out=-tt*in[1];
        }
        inline double sol0(double tt) {return cos(0.5*tt*tt);}
        inline double sol1(double tt) {return sin(0.5*tt*tt);}
};

class osc_rk4d : public rk4d, public osc {
    public:
        osc_rk4d() : rk4d(2) {}
        ~osc_rk4d() {}
        virtual void init() {osc_init(q);}
        virtual void ff(double tt,double *in,double *out) {osc_ff(tt,in,out);}
};

class osc_fsal : public fsal, public osc {
    public:
        osc_fsal() : fsal(2) {}
        ~osc_fsal() {}
        virtual void init() {osc_init(q);}
        virtual void ff(double tt,double *in,double *out) {osc_ff(tt,in,out);}
};

#endif