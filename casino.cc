#include <cmath>
#include <cstdio>
#include "custom_rng.hh"
#include "omp.h"

double play(custom_rng c){
	double x=0.;
	double win=0.;
	while(x<=1.){
		x+=c->doub();
		win+=100.;
	}
	return win;
}

int main(){

}