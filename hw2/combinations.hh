#ifndef COMBINATIONS_HH
#define COMBINATIONS_HH

class combinations{
	public:
		int count;
		int nn;
		int* all_comb;

		combinations(int nn_);
		~combinations();

		void helper(int* x,int r,int a, int b,int idx);
		int find_combination(int n,int r);
		void reset();


};

#endif