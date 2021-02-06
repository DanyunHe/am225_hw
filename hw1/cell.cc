#include <cmath>
#include <cstdio>
#include <cstring>

const int m=80;
const int n=40;
int* c=new int[m*n]; //cellular automaton on m by n grid 
int* cc=new int[m*n]; // next step cell grid 
custom_rng rng=custom_rng(0); // random number generator

// Initialize cell grid at time 0
void init(){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++,c++){
			if((i>=m/2-6)&&(i<m/2-6)&&(j>=n/2-6)&&(j<n/2-6)){
				if(rng->doub()>=3/4.){*c=1;}
				else{*c=0;}
			}
			else{*c=0;}
		}
	}
}

// Perform a gerneation step
void generate(){

	int nij;
	for(int i=0;i<m;i++){
		for(int j=0,c=0,cc=0;j<n;j++,c++,cc++){

			// Calculate Nij for cell at (i,j)
			nij=0;
			if(*(c+1)==1) nij++;
			if(*(c-1)==1) nij++;
			if(*(c+n)==1) nij++;
			if(*(c-n)==1) nij++;
			if(*(c+n+1)==1) nij++;
			if(*(c+n-1)==1) nij++;
			if(*(c-n+1)==1) nij++;
			if(*(c-n-1)==1) nij++;

			// Update cell status 
			if(*c==1){
				if(nij>=6){*cc=0;}
				else{*cc=1;}
			}
			else if(nij==3){*cc=1;}
			else {*cc=0;}
		}
	}

	memcpy(c[0],cc[0],m*n*sizeof(int));
}

// Save the cell information 
void output(char[] filename){

	FILE *fp=fopen("filename","wb");
	if(fp==NULL){
		fputs("error: fail to open file",stderr);
		exit(-1);
	}

	fwrite(c,sizeof(int),n*m,fp);
	fclose(fp);

}

int main(){

	init();
	for(int t=0;t<150;t++){
		if(t%25==0){
			char fn[10];
			sprintf(fn,"%d.out",t);
			output(fn);
		}
		gerenate();
	}

	delete [] c;
	delete [] cc;
	return 0;

}