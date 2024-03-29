#ifndef RBF_HH
#define RBF_HH

#include <cstdlib>

#include "conj_grad.hh"

// Data structure of (x,y) positions, with an additional variable z for sorting
// on
struct xyz {
    /** The x coordinate. */
    double x;
    /** The y coordinate. */
    double y;
    /** An additional value to sort on. */
    int z;
    xyz() {}
    xyz(double x_,double y_,int z_) : x(x_), y(y_), z(z_) {}
    xyz(const xyz &o) : x(o.x), y(o.y), z(o.z) {}
    void set(double x_,double y_,int z_) {
        x=x_;y=y_;z=z_;
    }
    const bool operator<(const xyz &o) {
        return z<o.z;
    }
    const bool operator==(const xyz &o) {
        return z==o.z;
    }
};

class rbf : public conj_grad {
    public:
        /** Number of specified points. */
        const int n;
        /** The type of radial function to use. */
        const int type;
        /** The square length scale of the radial function. */
        double lsq;
        /** The inverse square length scale of the radial function. */
        double ilsq;
        /** x positions of points. */
        double* const px;
        /** y positions of points. */
        double* const py;
        /** Function values at points. */
        double* const pf;
        /** Weights of the points in the RBF interpolant. */
        double* const rs;
        xyz* q;

        rbf(int n_,int type_);
        ~rbf();
        void init_random(int mode);
        void eigenvalues();
        void solve_weights_lapack();
        int solve_weights_conj_grad(int bls=0,bool verbose=false);
        void output_points(const char* filename);
        void output_interpolant(const char* filename,int q,double L);
        virtual void mul_A(double *in,double *out);
        virtual void M_inv(double *in,double *out);
        void make_table();
        int count_tk();
        int count_pk();
        inline void set_length_scale(double lscale) {
            lsq=lscale*lscale;ilsq=1./lsq;
        }
        inline int assemble_matrix() {
            if(A==NULL){
                int count=fill_matrix_entries(A=new double[n*n],0,n);
                return count;
            }
            return 0;
        }
    protected:
        double urand();
    private:
        int preconditioning_table(int bls_);
        int fill_matrix_entries(double *Ap,int k,int b);
        double phi(double rsq);
        int hilbert_curve(double px,double py);
        int to_int(double x,int fac);
        void rot(int n,int *x,int *y, int rx, int ry);
        /** The block size in the block Jacobi preconditioner. */
        int bls;
        /** The last block size in the block Jacobi preconditioner. */
        int lbls;
        /** A pointer to the dense matrix entries. */
        double *A;
        /** A pointer to the block Jacobi preconditioning information. */
        double *Apre;
        /** A pointer to the pivoting information in the block Jacobi
         * preconditioner. */
        int *ipiv;
        /** The number of entries in each row of the sparse matrix. */
        int* An;
        /** The indices of non-zero entries in the sparse matrix. */
        int* phi_index;
        /** The values of non-zero entries in the sparse matrix. */
        double* phi_val;
};

#endif