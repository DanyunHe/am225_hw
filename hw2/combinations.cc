#include "combinations.hh"
#include <cstring>
#include <cmath>
#include <cstdio>


combinations::combinations(int nn_): count(0),all_comb(new int[nn_]){}

combinations::~combinations(){
	delete [] all_comb;
}

void combinations::helper(int* x,int r,int a, int b,int idx){
	// printf("a %d b %d idx %d r %d x %d\n",a,b,idx,r,x[0]);
	if(idx==r){
		count++;
		// puts("c");
		for(int i=0;i<r;i++){
			// Copy x to the set of combinations
			all_comb[i+r*count]=x[i];
			// printf("i %d x[i] %d size %d\n",i,x[i],count);
		}
		return;
	}
	else if(a<=b){
		// puts("b");
		x[idx]=a;
		helper(x,r,a+1,b,idx+1);
		// printf("a %d b %d idx %d x %d\n",a,b,idx,x[0]);
		helper(x,r,a+1,b,idx);

	}
	return;
}

// Return number of combinations 
// Find Cn^r
int combinations::find_combination(int n, int r){

	int* x=new int[r];
	// puts("a");
	helper(x,r,0,n-1,0);
	return count;
	// for(auto it=*combinations.begin();it!=*combinations.end();it++){
	// 	std::cout<<*it<<"\n";
	// }
	// cout<<combinations;
	
	// puts("d");
	// return combinations.size();


}