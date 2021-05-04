#ifndef COLORS_HH
#define COLORS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct color {
    /** The red channel. */
    double r;
    /** The green channel. */
    double g;
    /** The blue channel. */
    double b;
    
    inline void prd_bc(color &c) {
        r=c.r;g=c.g;b=c.b;
    }
    inline void no_slip(color &c) {
        r=-c.r;g=-c.g;b=-c.b;
    }
    /** Computes the maximum allowable timestep based on the CFL restriction
     * from the velocity stored in this class.
     * \param[in] (xsp,ysp) the horizontal and vertical grid spacings.
     * \return The reciprocal of the maximum allowable timestep. */
    inline double cfl(double xsp,double ysp,double zsp) {
        double rc=fabs(r)*xsp,gc=fabs(g)*ysp,bc=fabs(b)*zsp;
        return rc>gc?rc:gc;
    }
};

#endif
