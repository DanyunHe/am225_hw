#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <iostream>
#include "combinations.hh"
using namespace std;


bool isvalid(int* nt,int k){
	for(int i=0;i<k-1;i++){
		if(nt[i]>nt[i+1]){
			return false;
		}
	}
	return true;
}

/** Count all trees of order p. */
int count_tree(int p){
	if(p==1){
		return 1;
	}
	if(p==2){
		return 1;
	}
	int count=0;
	// k: number of edges from the root
	// when k=1
	count+=count_tree(p-1);
	int num_com;
	combinations comb(1000);
	// combinations comb(1000);

	// Consider edges number from 2 to p-1
	for(int k=2;k<p;k++){
		// puts("1");
		// printf("n %d r %d\n",p-2,k-1);
		// combinations comb(1000);
		comb.reset();
		num_com=comb.find_combination(p-2,k-1);
		// printf("num %d\n",num_com);
		// puts("2");
		// for(int i=0;i<num_com*(k-1);i++) printf("i %d com %d\n",i,comb.all_comb[i]);

		int all_n[k];

		for(int i=0;i<num_com;i++){
			int nt=comb.all_comb[i*(k-1)]+1;
			all_n[0]=nt;
			// int total_tree=count_tree(nt);
			// printf("i %d nt %d\n",i,nt);
			for(int j=1;j<k-1;j++){
				nt=comb.all_comb[j+i*(k-1)]-comb.all_comb[j-1+i*(k-1)];
				all_n[j]=nt;
				// printf("i %d nt %d\n",i,nt);
				// total_tree*=count_tree(nt);
				// count+=count_tree(a[j+i*num_com]);
			}
			nt=p-1-comb.all_comb[k-2+i*(k-1)]-1;
			all_n[k-1]=nt;
			int total_tree=0;
			if(isvalid(all_n,k)){
				total_tree=1;
				for(int i=0;i<k;i++){
					total_tree*=count_tree(all_n[i]);
				}
			}
			// printf("i %d nt %d\n",i,nt);
			// total_tree*=count_tree(nt);
			count+=total_tree;

		}

	}
	return count;

}

int main(){
	for(int p=1;p<10;p++){
		int sol=count_tree(p);
		printf("order %d has trees %d \n",p,sol);

	}
	

}





