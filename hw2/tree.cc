#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <iostream>
#include "combinations.hh"
using namespace std;

// Make sure the numbers in the array nt is increasing 
bool isvalid(int* nt,int k){
	for(int i=0;i<k-1;i++){
		if(nt[i]>nt[i+1]){
			return false;
		}
	}
	return true;
}

// Compute factorial function n!
int factorial(int n){
	if(n==1) return 1;
	return n*factorial(n-1);
}

// Compute combination function Cn^r
int Cnr(int n, int r){
	if(r>n){
		printf("invalid n %d and r %d\n",n,r);
		return -1;
	}
	else if(r==n){
		return 1;
	}
	if(n==1) return 1;
	if(r==1) return n;

	return Cnr(n-1,r)+Cnr(n-1,r-1);

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
	combinations comb(10000000);
	// combinations comb(1000);

	// Consider edges number from 2 to p-1
	for(int k=2;k<p;k++){
		// puts("1");
		// if(k==p-1) printf("n %d r %d\n",p-2,k-1);
		// combinations comb(1000);
		comb.reset();
		num_com=comb.find_combination(p-2,k-1);
		// if(p==15) printf("edge %d num %d\n",k,num_com);
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
				if(p==7){
					// printf("k %d p %d\n",k,p);
					for(int id=0;id<k;id++){
						// printf("n[%d] %d\n",id,all_n[id]);
					}

				}
				total_tree=1;

				// Check if there is identical ni, then perform combination with repitation
				int ri=1,ni=all_n[0],num_tree;
				for(int id=0;id<k-1;id++){
					if(all_n[id]<3||all_n[id+1]!=ni){
						if(ri>1){
							// puts("a");
							num_tree=count_tree(ni);
							//perform combination with repitation
							// total_tree*=(factorial(ri+num_tree-1)/(factorial(ri)*factorial(num_tree-1)));
							total_tree*=Cnr(ri+num_tree-1,ri);
							ri=1;
						}
						ni=all_n[id+1];
						total_tree*=count_tree(all_n[id]);
						
					}
					else{
						ri++;
						ni=all_n[id+1];
						// printf("ri %d ni %d\n",ri,ni);
					}
				}
				if(ri>1){
					// puts("b");
					num_tree=count_tree(ni);
					// total_tree*=(factorial(ri+num_tree-1)/(factorial(ri)*factorial(num_tree-1)));
					total_tree*=Cnr(ri+num_tree-1,ri);
				}
				else{
					total_tree*=count_tree(all_n[k-1]);
				}
				// if(p==7) printf("k %d total_tree %d\n",k,total_tree);
			}
			// if(p==15) printf("count %d total_tree %d\n",count,total_tree);
			// total_tree*=count_tree(nt);
			count+=total_tree;


		}

	}
	return count;

}

int main(){
	for(int p=1;p<16;p++){
		int sol=count_tree(p);
		printf("order %d has trees %d \n",p,sol);

	}

	// int result=Cnr(10,3);
	// printf("c(10,3) %d",result);
	

}





