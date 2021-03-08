#include <cmath>
#include <cstring>
#include <list>

/** Count all trees of order p. */
int count_tree(int p){
	if(p==1){
		return 1;
	}
	int count=0;
	// k: number of edges from the root
	// when k=1
	count+=count_tree(p-1);
	int num_com;
	for(int k=2;k<p-1;k++){
		int* a; // the order of tree at the kth edge
		num_com=find_combination(k,a);
		for(int i=0;i<num_com;i++){
			// update a[0-k]
			for(int j=0;j<k;j++){
				count+=count_tree(a[j+i*num_com]);
			}
		}

	}

}

// Return number of combinations 
// Find Cn^r
int find_combination(int n, int r, list<int> combinations){
	list<int> combinations;
	int* x=new int[r];
	helper(combinations,x,0,n-1,0);

}

void helper(list<int> combinations, int* x,int r,int a, int b,int idx){
	if(idx==r){
		for(int i=0;i<r;i++){
			// Copy x to the set of combinations
			combinations.push_back(x[i]);
		}
	}
	else if(a<=b){
		x[idx]=a;
		helper(combinations,x,r-1,a+1,b,idx+1);
		helper(combinations,x,r,a+1,b,idx);

	}
}