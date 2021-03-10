#ifndef DISC_HH
#define DISC_HH

#include "sol_rk4d.hh"
#include "fsal.hh"
#include <cmath>
/** An ODE system with discountinuity. */
class disc {
    public:
        inline void disc_init(double *q) {
            *q=1.;q[1]=0.;
        }
        inline void disc_ff(double tt,double *in,double *out) {
            double ax=(in[0]>=0)?in[0]:-in[0];
            double ay=(in[1]>=0)?in[1]:-in[1];
            if(ax>=ay){
                out[1]=*in;
                *out=0.;
            }
            else{
                out[1]=0.;
                *out=-in[1];
            }
        }
        inline double sol0(double tt) {
            // solution for t=48+e^-1
            return 1.;
        }
        inline double sol1(double tt) {
            return exp(-1.);
        }
};

class disc_rk4d : public rk4d, public disc {
    public:
        disc_rk4d() : rk4d(2) {}
        ~disc_rk4d() {}
        virtual void init() {disc_init(q);}
        virtual void ff(double tt,double *in,double *out) {disc_ff(tt,in,out);}
};

class disc_fsal : public fsal, public disc {
    public:
        disc_fsal() : fsal(2) {}
        ~disc_fsal() {}
        virtual void init() {disc_init(q);}
        virtual void ff(double tt,double *in,double *out) {disc_ff(tt,in,out);}
};

#endif