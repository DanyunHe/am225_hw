#include <cmath>
#include <cstdio>
#include <cstring>


/** Custom random number generator based of "Ran" routine in Numerical Recipes
 * by Press et al. */
class custom_rng{
    public:
        unsigned long a,b,c; 
        custom_rng(unsigned long seed): b(4101842887655102017L), c(1) {
            if(sizeof(unsigned long)<8) {
                fputs("Error: 'unsigned long' type too short\n",stderr);
                exit(1);
            }
            a=seed^b;int64();
            b=a;int64();
            c=b;int64();
        }

        unsigned long int64() {
            a=a*2862933555777941757L+7046029254386353087L;
            b^=b>>17;b^=b<<31;b^=b>>8;
            c=4294957665U*(c&0xffffffff)+(c>>32);
            unsigned long d=a^(a<<21);
            d^=d>>35;d^=d<<4;
            return (d+b)^c;
        }
        inline double doub(){
            return 5.42101086242752217E-20*int64();
        }
};

const int m=80;
const int n=40;
int* c=new int[m*n]; //cellular automaton on m by n grid 
int* cc=new int[m*n]; // next step cell grid 
custom_rng rng=custom_rng(0); // random number generator

// Initialize cell grid at time 0
void init(){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			if((i>=m/2-6)&&(i<m/2-6)&&(j>=n/2-6)&&(j<n/2-6)){
				if(rng.doub()>=3/4.){c[j+n*i]=1;}
				else{c[j+n*i]=0;}
			}
			else{c[j+n*i]=0;}
		}
	}
}

// Perform a gerneation step
void generate(){

	int nij;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){

			// Calculate Nij for cell at (i,j)
			nij=0;
			int idx=j+i*n;
			if(c[idx+1]==1) nij++;
			if(c[idx-1]==1) nij++;
			if(c[idx+n]==1) nij++;
			if(c[idx-n]==1) nij++;
			if(c[idx+n+1]==1) nij++;
			if(c[idx+n-1]==1) nij++;
			if(c[idx-n+1]==1) nij++;
			if(c[idx-n-1]==1) nij++;

			// Update cell status 
			if(c[idx]==1){
				if(nij>=6){cc[idx]=0;}
				else{cc[idx]=1;}
			}
			else if(nij==3){cc[idx]=1;}
			else {cc[idx]=0;}
		}
	}

	memcpy(c,cc,m*n*sizeof(int));
}

// Save the cell information 
void output(char* filename){

	FILE *fp=fopen(filename,"wb");
	if(fp==NULL){
		fputs("error: fail to open file",stderr);
		exit(-1);
	}

	fwrite(c,sizeof(int),n*m,fp);
	fclose(fp);

}

int main(){

	init();
	for(int t=0;t<151;t++){
		if(t%25==0){
			char fn[10];
			sprintf(fn,"%d.out",t);
			output(fn);
		}
		generate();
	}

	delete [] c;
	delete [] cc;
	return 0;

}