#ifndef BILINEAR_2D_FE_HH
#define BILINEAR_2D_FE_HH

#include "conj_grad.hh"

#include <cstdio>
#include <cmath>

/** Class for solving an elliptic PDE problem over the domain [1,2] using
 * piecewise cubic basis functions. */
class bilinear_2d_fe: public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n;
        const int dof;
        /** The source function. */
        double* const f;
        /** The grid spacing. */
        const double h;
        /** The Neumann condition to apply at x=2. */
        double g;
        bilinear_2d_fe(int n_) : conj_grad((n_-1)*(n_-1)),
           n(n_),dof(n*n), f(new double[(n+1)*(n+1)]), h(2./n) {}
        virtual ~bilinear_2d_fe() {delete [] f;}
        void init();
        double l2_norm_mms();
        void print_matrix();
        void print(FILE *fp);
        void print(const char* filename);
        inline void print() {print(stdout);}
        virtual void mul_A(double *in,double *out);
    private:
        inline double mms_dsq(double vv,double ww,double s) {
            double del=exp(1-vv)*(3.+(vv-4.)*vv-ww*ww)-s;
            return del*del;
        }
        double quadracture_calc(double x2,double y2,double j2,double k2,double x3,double y3,double j3,double k3);
        void assemble_b();
        void cal_D(double xx,double yy,double* D);
        int basis_func(int j,int k);
        double basis_func_val(int idx,double xx,double yy);
        void dx_p(int j, int k,double xx, double yy,double &dx,double &dy);
        double integrand(double xx, double yy,int k,int j,int j2,int k2,int j3,int k3);
        double quadracture_calc(int k,int j,int j2,int k2,int j3,int k3);
        double quadracture_calc2(int k,int j,int j2,int k2,int j3,int k3);
        double integrand2(double xx, double yy,int k,int j,int j2,int k2,int j3,int k3);

};

#endif
