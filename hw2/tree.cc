#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <iostream>
#include "combinations.hh"
using namespace std;



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
	// combinations comb(1000);
	for(int k=2;k<p-1;k++){
		puts("1");
		printf("n %d r %d\n",p-1,k-1);
		combinations comb(1000);
		num_com=comb.find_combination(p-2,k-1);
		printf("num %d\n",num_com);
		puts("2");
		for(int i=0;i<num_com;i++){
			int nt=comb.all_comb[i*(k-1)]+1;
			int total_tree=count_tree(nt);
			for(int j=1;j<k-1;j++){
				nt=comb.all_comb[j+i*(k-1)]-comb.all_comb[j+1+i*(k-1)];
				total_tree*=count_tree(nt);
				// count+=count_tree(a[j+i*num_com]);
			}
			nt=p-1-comb.all_comb[k-2+i*(k-1)];
			total_tree*=count_tree(nt);

		}

	}
	return count;

}

int main(){
	int p=7;
	int sol=count_tree(p);
	printf("order %d has trees %d \n",p,sol);


}





