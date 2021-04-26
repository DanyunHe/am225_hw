#include <cstdlib>

#include "bilinear_2d_fe.hh"
#include "blas.h"
#include "quadrat.hh"

/** Initializes the source function to be a constant. */
void bilinear_2d_fe::init_const() {
    for(int i=0;i<=n*n;i++) f[i]=1.;
    assemble_b();
}

/** Initializes the source function to be a linear slope, f(x)=1.5-x. */
void bilinear_2d_fe::init_slope() {
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=xx-1.5;
    assemble_b();
}

/** Initializes the source function so that the solution will match
 * a manufactured solution, u(x)=exp(1-x)*sin(5*pi*x). */
void bilinear_2d_fe::init_mms() {
    const double o=5*M_PI;
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=-exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx));
    assemble_b();
}

/** Prints the stiffness matrix as a text array. */
void bilinear_2d_fe::print_matrix() {
    int i,j;

    // Allocate zero vector, and workspace to compute a matrix product
    double *r=new double[6*n],*s=r+3*n;
    for(i=0;i<3*n;i++) r[i]=0.;

    for(int i=0;i<3*n;i++) {

        // Apply the black box matrix multiplication routine to a unit vector,
        // in order to extract a column of matrix entries
        r[i]=1.;mul_A(r,s);r[i]=0.;

        // Print a row of matrix entries. This assumes the matrix is symmetric
        // (as required for conjugate gradient) so that the row<->column switch
        // is permissible.
        for(j=0;j<3*n-1;j++) printf("%g ",s[j]);
        printf("%g\n",s[3*n-1]);
    }
    delete [] r;
}

/** Performs multiplication on a vector by the stiffness matrix. */
void bilinear_2d_fe::mul_A(double *in,double *out) {
    int i,j,k;

    // Pre-computed integrals of derivatives of Lagrange polynomials, which are
    // required to construct the stiffness matrix
    const double B[16]={37/30.,-63/40.,9/20.,-13/120.,
                        -63/40.,18/5.,-99/40.,9/20.,
                        9/20.,-99/40.,18/5.,-63/40.,
                        -13/120.,9/20.,-63/40.,37/30.},
                 C[16]={17/40.,-51/80.,3/8.,-13/80.,
                        -51/80.,27/8.,-297/80.,39/40.,
                        3/8.,-297/80.,297/40.,-327/80.,
                        -13/80.,39/40.,-327/80.,131/40.};

    // Set the output vector to initially be zero
    for(i=0;i<dof;i++) out[i]=0.;

    // Loop over the square elements
    for(k=0;k<n;k++) for(int j=0;j<n;j++) {

        int jlo=j==0?1:0,jhi=j==n-1?0:1,
            klo=k==0?1:0,kji==hi=k==n-1?0:1,
            j2,k2,j3,k3;

        // Loop over different basis functions in this interval. Here i is
        // alpha, and j is beta as defined in the notes. Usually i and j run
        // from 0 to 3. For the first interval when k=0, the first basis
        // function is not included because it is not present due to the
        // essential boundary condition. Hence i and j run from 1 to 3 in that
        // case.
        double result;
        for(k2=0;k2<=1;k2++) for(j2=0;j2<=1;j2++){
            for(k3=klo;k3<=khi;k3++) for(j3=jlo;j3<=jhi;j3++){
                
                // Compute contributaion for (j,k) square
                // Using 1st bilinear func indexed with (j2,k2);
                // Using 2nd bilinear func indexed with (j3,k3);
                result=quadracture_calc(x2,y2,j2,k2,x3,y3,j3,k3);

                // Store into location (j+j2-1,k+k2-1)
                ind2=(j+j2-1)+(n-1)*(k+k2-1);
                ind3=(j+j3-1)+(n-1)*(k+k3-1);
                out[ind3]+=result*in[ind2];

            }

        }
    }
}

void cal_D(double x,double y,double* D){
    D[0]=sqrt(1.-y*y/2);
    D[1]=-x*y/(2.*sqrt(1.-y*y/2.));
    D[2]=-x*y/(2.*sqrt(1.-x*x/2.));
    D[3]=sqrt(1.-x*x/2);
}

void dx_p(int j, int k,double x, double y,double &dx,double &dy){
    if(j==0){
        if(k==0){
            dx=-1+y;
            dy=-1+x;
        }
        else if(k==1){
            dx=-y;
            dy=1-x;
        }
    }
    else if(k==0){
        dx=1-y;
        dy=-x;

    }
    else{
        dx=y;
        dy=x;
    }

}

/** x,y in range [0,1]. */
double f(double x, double y,int k,int j,int j2,int k2,int j3,int k3){
    double* D=new double[4];
    cal_D(-1+h*(x+k),-1+h*(y+l),D);
    double dx2,dy2,dx3,dy3;
    dx_p(j2,k2,x,y,dx2,dy2);
    dx_p(j3,k3,x,y,dx3,dy3);
    result=inv(D)*dx2*dy2*D*det(D);
    return result;

}

double quadracture_calc(int k,int j,int j2,int k2,int x3,int y3,int j3,int k3){
    // Set up the quadrature points and weights
    int np=5;
    double I,Ir;
    quadrat q(np);

    // Perform the 2D sum of function evaluations, each multiplied by the
    // corresponding weight
    I=0.;
    for(int b=0;b<np;b++) {
        Ir=0.;
        for(int a=0;a<np;a++) Ir+=q.w[a]*f(q.x[a],q.x[b],k,j,j2,k2,j3,k3);
        I+=Ir*q.w[b];
    }
    return I;

}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void bilinear_2d_fe::assemble_b() {
    int i,j,k;

    // Pre-computed integrals of Lagrange polynomial products
    const double D[16]={8/35.,99/560.,-9/140.,19/560.,
                        99/560.,81/70.,-81/560.,-9/140.,
                        -9/140.,-81/560.,81/70.,99/560.,
                        19/560.,-9/140.,99/560.,8/35.};

    // Clear the source function
    for(i=0;i<3*n;i++) b[i]=0.;

    // Loop over each interval. Here k=3q as defined in the notes.
    for(k=0;k<3*n;k+=3) {

        // Loop over different basis functions in this interval. Here
        // i is alpha and j is beta as defined in the notes. Usually i and j
        // run from 0 to 3, but for the first interval i only runs from 1 to 3
        // because of the essential boundary condition.
        for(i=(k==0?1:0);i<4;i++) for(j=0;j<4;j++) {

            // Add a contribution to the source vector from the two basis
            // functions. The -1 shift is included when indexing b since the
            // solution at x_0=1 is not part of the vector. However, f does
            // contain a term at x_0=1, and thus there is no shift applied
            // there.
            b[-1+k+i]+=D[i+4*j]*f[k+j];
        }
    }

    // Normalize the results, and add in the Neumann condition to the last
    // entry
    for(i=0;i<3*n;i++) b[i]*=h;
    b[3*n-1]+=2*g;
}

/** Prints the solution.
 * \param[in] fp a file handle to write to. */
void bilinear_2d_fe::print(FILE *fp) {
    double xx=1+h;
    fprintf(fp,"1 0 0 %g\n",*f);
    for(int i=0;i<3*n;i++,xx+=h) fprintf(fp,"%g %g %g %g\n",xx,x[i],b[i],f[i+1]);
}

/** Prints the solution.
 * \param[in] filename the name of the file to write to. */
void bilinear_2d_fe::print(const char* filename) {
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
double bilinear_2d_fe::l2_norm_mms() {
    double l2=0.,xx=1.;
    for(int i=0;i<3*n-1;i++) {
        xx+=h;
        l2+=mms_dsq(xx,x[i]);
    }

    // Add the contribution at the last point, including a factor of 1/2 due to
    // the trapezoid rule. Note that there is no contribution at x=1, since the
    // numerical solution is zero there, and hence matches the manufactured
    // solution perfectly.
    l2+=0.5*mms_dsq(2.,x[3*n-1]);
    return sqrt(h*l2);
}