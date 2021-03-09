#ifndef ORBIT_HH
#define ORBIT_HH

#include "sol_geng.hh"

#include <cstdio>
#include <cmath>

/** This class has functions to specify the test orbit problem. */
class orbit {
    public:
        /** The initial value of the Hamiltonian, which is used to track
         * changes in the value over time. */
        double init_h;
        /** Computes the change in the p (momentum) variables from the q
         * (coordinate) variables.
         * \param[in] in the current state of the p and q variables.
         * \param[in] out the change in the p variables. */
        inline void orb_fq(double *in,double *out) {
            const double Omega = 0.25;
            double q1 = in[3];
            double q2 = in[4];
            out[3]=in[0]+q2*Omega;
            out[4]=in[1]-q1*Omega;
            out[5]=in[2];
        }
        /** Computes the change in the q (momentum) variables from the p
         * (coordinate) variables.
         * \param[in] in the current state of the p variables.
         * \param[in] out the change in the q variables. */
        inline void orb_fp(double *in,double *out) {
            const double Omega = 0.25;
            const double A = 1.0;
            const double C = 1.0;
            const double a = 1.25;
            const double b = 1.0;
            const double c = 0.75;
            double B = C+(in[3]*in[3])/(a*a)+(in[4]*in[4])/(b*b)+(in[5]*in[5])/(c*c);
            double p1 = in[0];
            double p2 = in[1];
            out[0]=Omega*p2-(2*A/(B*a*a))*in[3];
            out[1]=-Omega*p1-(2*A/(B*b*b))*in[4];
            out[2]=-(2*A/(B*c*c))*in[5];
        }
        void orb_init(double *q);
        void orb_print(double t_,double *q);
        double hamiltonian(double *q);
};

/* Class to solve using Geng's method*/
class orb_geng : public geng, public orbit {
    public:
        orb_geng(double a_, double e_) : geng(6), orbit() {}
        virtual void ff(double t_,double *in,double *out) {
            orb_fp(in,out);orb_fq(in,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

#endif
