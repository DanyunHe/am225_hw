#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "fluid_2d.hh"

const char fn[]="ftest.out";

int main() {

    // Create the output directory for storing the simulation frames
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Specify which fields should be outputted. 1: horizontal velocity, 2:
    // vertical velocity, 4: pressure.
    unsigned int fflags=1|2|4;

    // Construct the simulation class, setting the number of gridpoints, the
    // periodicity, and physical constants
    fluid_2d f2d(1024,1024,false,false,-1,1,-1,1,0.002,1.,fflags,fn);

    // Initialize the tracers, and set the timestep based on multiplying the
    // maximum allowable by a padding factor
    f2d.initialize(8,0.3);

    // Run the simulation for a specified duration, outputting snapshots at
    // regular intervals
    f2d.solve(0.2,10);
    f2d.add_layer();
    f2d.solve(0.5,10);
}
