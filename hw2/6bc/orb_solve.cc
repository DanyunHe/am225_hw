#include <cstdio>
#include <cmath>
#include <cstring>

#include "sol_geng.hh"
#include "orbit.hh"

int main(int argc,char **argv) {

    // Number of timesteps
    const int steps=1e5*20;

    // The semi-major axis of the orbit
    const double a=1.;

    // The eccentricity of the orbit
    const double e=1/3.;

    // The simulation duration, currently set to three complete orbits
    const double du=1e5;//600*M_PI*a*sqrt(a);

    // Check for one command-line argument
    if(argc!=2) {
        fputs("Syntax: ./orb_solve <method_type>\n\n"
              "Solves an elliptical orbit problem with a variety of non-symplectic (NS) and\n"
              "symplectic (S) methods.\n\n"
              "Types:\n"
              "'geng' - Geng's method, S\n",stderr);
        return 1;
    }

    // Check the method type, and call the corresponding integration routine
    if(strcmp(argv[1],"geng")==0) {
    	// Geng
        orb_geng o(a,e);o.solve_fixed(du,steps,true);
    } else fputs("Invalid method type\n",stderr);
}
