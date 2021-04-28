#include <cstdlib>

#include "alter_1d_fe.hh"
#include "blas.h"

/** Initializes the source function to be a constant. */
void alter_1d_fe::init_const() {
    for(int i=0;i<=3*n;i++) f[i]=1.;
    assemble_b();
}

/** Initializes the source function to be a linear slope, f(x)=1.5-x. */
void alter_1d_fe::init_slope() {
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=xx-1.5;
    assemble_b();
}

/** Initializes the source function so that the solution will match
 * a manufactured solution, u(x)=exp(1-x)*sin(5*pi*x). */
void alter_1d_fe::init_mms() {
    const double o=5*M_PI;
    double xx=1;
    for(int i=0;i<=n;i++,xx+=h) {
        f[2*i]=-exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx));
        f[2*i+1]=exp(1-xx)*(o*(4+(-3+o*o)*xx)*cos(o*xx)+(-2+o*o*(2-3*xx)+xx)*sin(o*xx))*h;
    }
    assemble_b();
}

/** Prints the stiffness matrix as a text array. */
void alter_1d_fe::print_matrix() {
    int i,j;

    // Allocate zero vector, and workspace to compute a matrix product
    double *r=new double[4*n+2],*s=r+2*n+1;
    for(i=0;i<2*n+1;i++) r[i]=0.;

    for(int i=0;i<2*n+1;i++) {

        // Apply the black box matrix multiplication routine to a unit vector,
        // in order to extract a column of matrix entries
        r[i]=1.;mul_A(r,s);r[i]=0.;

        // Print a row of matrix entries. This assumes the matrix is symmetric
        // (as required for conjugate gradient) so that the row<->column switch
        // is permissible.
        for(j=0;j<2*n;j++) printf("%g ",s[j]);
        printf("%g\n",s[2*n]);
    }
    delete [] r;
}

/** Performs multiplication on a vector by the stiffness matrix. */
void alter_1d_fe::mul_A(double *in,double *out) {
    int i,j,k;

    // Pre-computed integrals of derivatives of Lagrange polynomials, which are
    // required to construct the stiffness matrix
    const double B[16]={1.2,0.1,-1.2,0.1,
                        0.1,2/15.,-0.1,-1/30.,
                        -1.2,-0.1,1.2,-0.1,
                        0.1,-1/30.,-0.1,2/15.},
                 C[16]={0.6,0.1,-0.6,0.,
                        0.1,1/30.,-0.1,-1/60.,
                        -0.6,-0.1,0.6,0.,
                        0.,-1/60.,0.,0.1};

    // Set the output vector to initially be zero
    for(i=0;i<2*n+1;i++) out[i]=0.;

    // Loop over each interval, and compute the contribution from
    // each
    for(k=0;k<2*n;k+=2) for(i=(k==0?1:0);i<4;i++)
        for(j=(k==0?1:0);j<4;j++)
            out[-1+k+i]+=((k/2.+1./h)*B[i+4*j]+C[i+4*j])*in[-1+k+j];
}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void alter_1d_fe::assemble_b() {
    int i,j,k;

    // Pre-computed integrals of Lagrange polynomial products
    const double D[16]={13/35.,11/210.,9/70.,-13/420.,
                        11/210.,1/105.,13/420.,-1/140.,
                        9/70.,13/420.,13/35.,-11/210.,
                        -13/420.,-1/140.,-11/210.,1/105.};

    // Clear the source function
    for(i=0;i<2*n+1;i++) b[i]=0.;

    // Loop over each interval, and compute the contributions
    // from each Lagrange polynomial pair
    for(k=0;k<2*n;k+=2) for(i=(k==0?1:0);i<4;i++)
        for(j=0;j<4;j++) b[-1+k+i]+=D[i+4*j]*f[k+j];

    // Normalize the results, and add in the Neumann condition to the last
    // entry
    for(i=0;i<2*n+1;i++) b[i]*=h;
    b[2*n-1]+=2*g;
}

/** Prints the solution.
 * \param[in] fp a file handle to write to. */
void alter_1d_fe::print(FILE *fp) {
    double xx=1+h;
    fprintf(fp,"1 0 0 %g\n",*f);
    for(int i=0;i<n;i++,xx+=h){
        fprintf(fp,"%g %g %g %g\n",xx,x[2*i],b[2*i],f[2*i+1]);
        fprintf(fp,"%g %g %g %g\n",xx,x[2*i+1],b[2*i+1],f[2*i+2]);
    }
}

/** Prints the solution.
 * \param[in] filename the name of the file to write to. */
void alter_1d_fe::print(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    print(fp);
    fclose(fp);
}

/** Computes the L2 norm between the numerical solution and the true
 * manufactured solution. The routine uses the trapezoid rule to evaluate the
 * integral.
 * \return The L2 norm. */
double alter_1d_fe::l2_norm_mms() {
    double l2=0.,xx=1.;
    for(int i=0;i<n;i++) {
        xx+=h;
        l2+=mms_dsq(xx,x[2*i+1]);
    }

    // Add the contribution at the last point, including a factor of 1/2 due to
    // the trapezoid rule. Note that there is no contribution at x=1, since the
    // numerical solution is zero there, and hence matches the manufactured
    // solution perfectly.
    l2+=0.5*mms_dsq(2.,x[2*n+1]);
    return sqrt(h*l2);
}
