#include <cstdlib>
#include <cmath>

#include "bilinear_2d_fe.hh"
#include "blas.h"
#include "quadrat.hh"

/** Initializes the source function to be a constant. */
void bilinear_2d_fe::init() {
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            double xx=-1.+h*i,yy=-1.+h*j;
            double v=xx*sqrt(1.-yy*yy/2.),w=yy*sqrt(1.-xx*xx/2.);
            f[j+(n+1)*i]=exp(-v)*(3.+(v-4.)*v+w*w);
        }
    } 
    assemble_b();
}


/** Prints the stiffness matrix as a text array. */
void bilinear_2d_fe::print_matrix() {
    int i,j;

    // Allocate zero vector, and workspace to compute a matrix product
    double *r=new double[2*dof],*s=r+dof;
    for(i=0;i<dof;i++) r[i]=0.;

    for(int i=0;i<dof;i++) {

        // Apply the black box matrix multiplication routine to a unit vector,
        // in order to extract a column of matrix entries
        r[i]=1.;mul_A(r,s);r[i]=0.;

        // Print a row of matrix entries. This assumes the matrix is symmetric
        // (as required for conjugate gradient) so that the row<->column switch
        // is permissible.
        for(j=0;j<dof-1;j++) printf("%g ",s[j]);
        printf("%g\n",s[dof-1]);
    }
    delete [] r;
}

/** Performs multiplication on a vector by the stiffness matrix. */
void bilinear_2d_fe::mul_A(double *in,double *out) {
    int i,j,k;

    // Set the output vector to initially be zero
    for(i=0;i<dof;i++) out[i]=0.;

    // Loop over the square elements
    for(k=0;k<n;k++) for(int j=0;j<n;j++) {

        int jlo=j==0?1:0,jhi=j==n-1?0:1,
            klo=k==0?1:0,khi=k==n-1?0:1,
            j2,k2,j3,k3;

        // Loop over different basis functions in this interval. Here i is
        // alpha, and j is beta as defined in the notes. Usually i and j run
        // from 0 to 3. For the first interval when k=0, the first basis
        // function is not included because it is not present due to the
        // essential boundary condition. Hence i and j run from 1 to 3 in that
        // case.
        double result;
        int ind2,ind3;
        for(k2=klo;k2<=khi;k2++) for(j2=jlo;j2<=jhi;j2++){
            for(k3=klo;k3<=khi;k3++) for(j3=jlo;j3<=jhi;j3++){
                
                // Compute contributaion for (j,k) square
                // Using 1st bilinear func indexed with (j2,k2);
                // Using 2nd bilinear func indexed with (j3,k3);
                result=quadracture_calc(j,k,j2,k2,j3,k3);
                // printf("result %g\n",result);

                // Store into location (j+j2-1,k+k2-1)
                ind2=(j+j2-1)+(n-1)*(k+k2-1);
                ind3=(j+j3-1)+(n-1)*(k+k3-1);
                // printf("ind3 %d ind2 %d\n",ind3,ind2);
                out[ind3]+=result*in[ind2];

            }

        }
    }
}

void bilinear_2d_fe::cal_D(double xx,double yy,double* D){
    D[0]=sqrt(1.-yy*yy/2);
    D[1]=-xx*yy/(2.*sqrt(1.-yy*yy/2.));
    D[2]=-xx*yy/(2.*sqrt(1.-xx*xx/2.));
    D[3]=sqrt(1.-xx*xx/2);
}

/** Return the basis function index. */
int bilinear_2d_fe::basis_func(int j,int k){
    if(j==0){
        if(k==0){
            return 0;
        }
        else if(k==1){
            return 2;
        }
    }
    else if(k==0){
        return 1;

    }
    else{
        return 3;
    }

}

/** Return basis function values at (x,y) for function idx.
 * f0(x,y)=(1-x)*(1-y) at (0,0)
 * f1(x,y)=x*(1-y) at (1,0)
 * f2(x,y)=(1-x)*y at (0,1)
 * f3(x,y)=x*y at (1,1) */
double bilinear_2d_fe::basis_func_val(int idx,double xx,double yy){
    switch(idx){
        case 0:
            return (1-xx)*(1-yy);
            break;
        case 1:
            return xx*(1-yy);
            break;
        case 2:
            return (1-xx)*yy;
            break;
        case 3:
            return xx*yy;
            break;
    }

}

/** Calculate partial derivative for phi_j,k at (x,y). */
void bilinear_2d_fe::dx_p(int j, int k,double xx, double yy,double &dx,double &dy){
    int func_idx=basis_func(j,k);
    switch(func_idx){
        case 0:
            dx=-1+yy;
            dy=-1+xx;
            break;
        case 1:
            dx=1-yy;
            dy=-xx;
            break;
        case 2:
            dx=-yy;
            dy=1-xx;
            break;
        case 3:
            dx=yy;
            dy=xx;
            break;
    }

}

/** Calculate integrand for x,y in range [0,1]. */
double bilinear_2d_fe::integrand(double xx, double yy,int j,int k,int j2,int k2,int j3,int k3){
    double* D=new double[4];
    cal_D(-1.+h*(xx+j),-1.+h*(yy+k),D);
    double dx2,dy2,dx3,dy3;
    dx_p(j2,k2,xx,yy,dx2,dy2);
    dx_p(j3,k3,xx,yy,dx3,dy3);

    double det_D=D[0]*D[3]-D[1]*D[2];
    // Calculate inverse of D
    double* inv_D=new double[4];
    inv_D[0]=D[3]/det_D,inv_D[1]=-D[1]/det_D,inv_D[2]=-D[2]/det_D,inv_D[3]=D[0]/det_D;
    double result=(inv_D[0]*dx2+inv_D[2]*dy2)*(inv_D[0]*dx3+inv_D[2]*dy3)+\
        (inv_D[1]*dx2+inv_D[3]*dy2)*(inv_D[1]*dx3+inv_D[3]*dy3);
    result*=det_D;
    return result;

}

double bilinear_2d_fe::quadracture_calc(int j,int k,int j2,int k2,int j3,int k3){
    // Set up the quadrature points and weights
    int np=5;
    double I,Ir;
    quadrat q(np);

    // Perform the 2D sum of function evaluations, each multiplied by the
    // corresponding weight

    // this from -1 to 1? adjust to 0 to 1
    I=0.;
    for(int jj=0;jj<np;jj++) {
        Ir=0.;
        for(int ii=0;ii<np;ii++){
            Ir+=q.w[ii]*0.5*integrand(q.x[ii]*0.5+0.5,q.x[jj]*0.5+0.5,j,k,j2,k2,j3,k3);
        }
        I+=Ir*q.w[jj]*0.5;
    }
    return I;

}



/** Calculate integral value for assemble_b function. */
double bilinear_2d_fe::quadracture_calc2(int k,int j,int j2,int k2,int j3,int k3){
    // Set up the quadrature points and weights
    int np=5;
    double I,Ir;
    quadrat q(np);

    // Perform the 2D sum of function evaluations, each multiplied by the
    // corresponding weight

    // this from -1 to 1? adjust to 0 to 1
    I=0.;
    for(int jj=0;jj<np;jj++) {
        Ir=0.;
        for(int ii=0;ii<np;ii++){
            Ir+=0.5*q.w[ii]*integrand2(q.x[ii]*0.5+0.5,q.x[jj]*0.5+0.5,k,j,j2,k2,j3,k3);
        }
        I+=Ir*q.w[jj]*0.5;
    }
    return I;

}

/** Calculate integrand for assemble_b function. */
double bilinear_2d_fe::integrand2(double xx, double yy,int k,int j,int j2,int k2,int j3,int k3){
    double* D=new double[4];
    cal_D(-1+h*(xx+k),-1+h*(yy+j),D);
    // Get basis function index
    int p2=basis_func(j2,k2);
    int p3=basis_func(j3,k3);

    double det_D=D[0]*D[3]-D[1]*D[2];
    double result=basis_func_val(p2,xx,yy)*basis_func_val(p3,xx,yy)*det_D;
    return result;

}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void bilinear_2d_fe::assemble_b() {
    int i,j,k;

    // Set the output vector to initially be zero
    for(i=0;i<dof;i++) b[i]=0.;

    // Loop over the square elements
    for(k=0;k<n;k++) for(int j=0;j<n;j++) {

        int jlo=j==0?1:0,jhi=j==n-1?0:1,
            klo=k==0?1:0,khi=k==n-1?0:1,
            j2,k2,j3,k3;

        // Loop over different basis functions in this interval. Here i is
        // alpha, and j is beta as defined in the notes. Usually i and j run
        // from 0 to 3. For the first interval when k=0, the first basis
        // function is not included because it is not present due to the
        // essential boundary condition. Hence i and j run from 1 to 3 in that
        // case.
        double result;
        int ind2,ind3;
        for(k2=klo;k2<=khi;k2++) for(j2=jlo;j2<=jhi;j2++){
            for(k3=0;k3<=1;k3++) for(j3=0;j3<=1;j3++){
                

                // Store into location (j+j2-1,k+k2-1)
                ind2=(j+j2-1)+(n-1)*(k+k2-1);
                ind3=j+j3+(n+1)*(k+k3);

                result=quadracture_calc2(k,j,j2,k2,j3,k3);
                b[ind2]+=result*f[ind3];
                // printf("ind2 %d result %g\n",ind2,result*f[ind3]);

            }

        }
    }
    // Normalize the results, and add in the Neumann condition to the last
    // entry
    for(i=0;i<dof;i++) b[i]*=(h*h);
    // for(i=0;i<dof;i++) printf("b[%d] %g\n",i,b[i]);
    //b[3*n-1]+=2*g;
}

/** Prints the solution.
 * \param[in] fp a file handle to write to. */
void bilinear_2d_fe::print(FILE *fp) {
    double xx,yy;
    fprintf(fp,"-1 -1 0 0 %g\n",*f);
    for(int i=0;i<n-1;i++){
        for(int j=0;j<n-1;j++){
            xx=-1.+(i+1)*h;
            yy=-1.+(j+1)*h;
            fprintf(fp,"%g %g %g %g %g\n",xx,yy,x[j+(n-1)*i],b[j+(n-1)*i],f[j+(n+1)*i]);

        }
    }
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
    double l2=0.,xx,yy,vv,ww;
    for(int j=0;j<n-1;j++) {
        for(int i=0;i<n-1;i++){
            xx=-1.+(i+1)*h;
            yy=-1.+(j+1)*h;
            vv=xx*sqrt(1.-yy*yy/2.);
            ww=yy*sqrt(1.-xx*xx/2.);
            l2+=mms_dsq(vv,ww,x[j+i*(n-1)]);
            // printf("xx %g yy %g l2 %g\n",xx,yy,l2);
        }
        // l2+=mms_dsq(vv,ww,x[j+j*n]);
    }

    // Add the contribution at the last point, including a factor of 1/2 due to
    // the trapezoid rule. Note that there is no contribution at x=1, since the
    // numerical solution is zero there, and hence matches the manufactured
    // solution perfectly.
    // l2+=0.5*mms_dsq(sqrt(1./2.),sqrt(1./2.),x[dof-1]);
    return sqrt(h*h*l2);
}
