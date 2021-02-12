#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include "omp.h"

// Check if x, the grid value at y, is a valid value to fill in board b
bool valid(int x,int y,int* b){

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
	int ic=ii-ii%3;
	int jc=jj-jj%3;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			int idx=(ic+i)*9+jc+j;
			if(b[idx]==x && idx!=y) return false;
		}
	}

	return true;

}


// Show the sudoku board
void print_board(int* b){
	for(int i=0;i<9;i++){
		for(int j=0;j<8;j++){
			printf("%d|",b[i*9+j]);
		}
		printf("%d\n",b[i*9+8]);
	}
}



unsigned long int SudokuSolve(int i,int* b){
	if(i==81){
		// puts("sudoku complete!");
		// print_board(b);
		return 1;
	}
	else if(b[i]!=0){
		return SudokuSolve(i+1,b);
	}
	else{
		unsigned long int total_sol=0;

	// Multithread searching 
	#pragma omp parallel
	{ 
		
		int* board=new int[81];
		memcpy(board,b,81*sizeof(int));

	#pragma omp for schedule(dynamic) reduction(+:total_sol)
		for(int j=1;j<10;j++){
			if(valid(j,i,board)){
				board[i]=j;
				total_sol+=SudokuSolve(i+1,board);
			}
		}
		board[i]=0;
		memcpy(b,board,81*sizeof(int));
		delete [] board;
		// return total_sol;
	}
	// if(total_sol>100000) printf("number of solutions for fig 1c %lu\n",total_sol);
	return total_sol;
	}
	
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
void board3(int* b){

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


	int b[81];
	int result;
	double t0,t1;
	// 4a
	int nrepeat=1;
	t0=omp_get_wtime();
	for(int i=0;i<nrepeat;i++){
		board1(b);
	// print_board(b1);
		result=SudokuSolve(0,b);

	}
	t1=omp_get_wtime();
	printf("time to solve a is %f\n",(t1-t0)/nrepeat);
	
	// // 4b
	t0=omp_get_wtime();
	board2(b);
	result=SudokuSolve(0,b);
	t1=omp_get_wtime();
	printf("number of solutions for fig 1b %d, time %f\n",result,t1-t0);

	// 4c
	t0=omp_get_wtime();
	board3(b);
	unsigned long int count=SudokuSolve(0,b);
	t1=omp_get_wtime();
	// std::cout<<"board c"<<std::dec<<count<<std::endl;
	printf("number of solutions for fig 1c %#lu, time %f\n",count,t1-t0);
	// printf("number of solutions for fig 1c %#lu, time %f\n",count,t1-t0);
	// printf("number of solutions for fig 1c %#lx, time %f\n",count,t1-t0);

}