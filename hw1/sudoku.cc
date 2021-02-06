#include <cmath>
#include <cstdio>

void SudokuSolve(int i,int* b){
	if(i==81){
		puts("sudoku complete!")
		for(int k=0;k<81;k++){
			printf("k %d, b[k] %d\n",k,b[k]);
		}
		return;
	}
	else{
		if(b[i]==0){
			for(int j=1;j<10;j++){
				if(valid(j,i,b)){
					b[i]=j;
					SudokuSolve(i+1);
					b[i]=0;
				}
			}
		}
		else{
			SudokuSolve(i+1);
		}
	}
}

// Check if x, the grid value at y, is a valid value to fill in board b
Boolean valid(int x,int y,int* b){

	// Get the row and column index (ii,jj)
	int ii=y/9;
	int jj=y%9;

	// Check row
	for(int j=0;j<9;j++){
		if(b[ii*9+j]==x) return false;
	}

	// Check column 
	for(int i=0;i<9;i++){
		if(b[i*9+jj]==x) return false;
	}

	// Check block
	// Get the index of the first grid in the block
	int ic=ii/3;
	int jc=jj/3;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			idx=(ic+i)*9+jc+j
			if(b[idx]==x && idx!=y) return false;
		}
	}

	return true;

}

// Initialize board in figure 1a
void board1(int* b){
	
	// indexing b[i,j]=b[i*9+j]
	for(int i=0;i<81;i++){b[i]=0;}
	
	b[0]=1,b[3]=9,b[5]=7,b[8]=3;
	b[10]=8,b[16]=7;
	b[20]=9,b[24]=6;

	b[29]=7,b[30]=2,b[32]=9,b[33]=4;
	b[36]=4,b[37]=1,b[43]=9,b[44]=5;
	b[47]=8,b[48]=5,b[50]=4,b[51]=3;

	b[56]=3,b[60]=7;
	b[64]=5,b[70]=4;
	b[72]=2,b[75]=8,b[77]=6,b[80]=9;
	
}

// Initialize board in figure 1b
void board2(int* b){

	// indexing b[i,j]=b[i*9+j]
	for(int i=0;i<81;i++){b[i]=0;}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			b[30*i+j]=j+1;
			b[30*i+9+j]=j+4;
			b[30*i+18+j]=j+7;
		}
	}
}

// Initialize board in figure 1c
void board2(int* b){

	// indexing b[i,j]=b[i*9+j]
	for(int i=0;i<81;i++){b[i]=0;}
	for(int i=0;i<3;i++){
		b[i*30]=1;
		b[i*30+10]=5;
		b[i*30+20]=9;

		b[i*24+8]=2;
		b[i*24+16]=4;
		b[i*24+24]=6;
	}
	b[40]=8;

}

int main(){

	int b1[81];
	board1(b1);
	SudokuSolve(0,b1);


}