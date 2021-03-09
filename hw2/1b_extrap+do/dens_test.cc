#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sol_rk4d.hh"
#include "osc.hh"
#include "brusselator.hh"

/** The number of components in the ODE system. */
const int ns=2;

int main() {

    brus_fsal o;
    o.solve_fixed(20.,2000,true,0,1e-10);
}
