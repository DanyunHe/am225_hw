#include <cmath>
#include <cstdio>
#include "custom_rng.hh"
#include "omp.h"

// Play a casino game and return the amount of winnings 
int play(custom_rng c){
	double x=0.;
	int win=0.;
	while(x<=1.){
		x+=c->doub();
		win+=100;
	}
	return win;
}

int main(){

	// Play 1e9 games and calculate the total winnings 
	int w=0;
#pragma omp parallel for reduction(+:w)
	for(int i=0;i<1e9;i++){
		custom_rng c=new custom_rng(i);
		w+=play(c);
		delete [] c;
	}

}