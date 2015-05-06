#include <stdlib.h>
#include <stdio.h>
#include "linear_solver.h"

int main() {

	srand (time(0));

	// Generate random Ax = b
	int n = 16;
	int row, col;
	double **A, *b, *v, *x;

	A = (double**) malloc (n*sizeof(double*));
	b = (double*) malloc (n*sizeof(double));
	v = (double*) malloc (n*sizeof(double));
	x = (double*) malloc (n*sizeof(double));
	
	// allocate and randomize
	for (row=0; row<n; row++) {
		A[row] = (double*) malloc(n*sizeof(double));
		for (col=0; col<n; col++) {
			A[row][col] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
		}
		b[row] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
	}

	// printf("prob. A="); print_sq_mat(n,A); printf("\n");
	// printf("prob. b="); print_vec(n,b); printf("\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	printf("conj. steps = %d\n", conjugate_solver (n, A, b, x));
	// validate result
	mat_vec_mult(n, A, x, v);
	printf("conj. error = %f\n",two_norm (n,b,v));
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	printf("stps. steps = %d\n", steepest_solver (n, A, b, x));
	// validate result
	mat_vec_mult(n, A, x, v);
	printf("stps. error = %f\n",two_norm (n,b,v));
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	for (row=0; row<n; row++) {
		free(A[row]);
	}
	free (A);
	free (b);
	free (v);
	free (x);

	return 0;

}