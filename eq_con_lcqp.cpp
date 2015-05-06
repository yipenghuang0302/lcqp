#include "lcqp.h"

double eq_con_lcqp (
	int n,
	double ** G,
	double * g0,
	int p,
	double ** CE, // p*n
	double * ce0, // p
	double * x,
	double * lagrange
) {

	int row, col;

	// allocate and construct K = [[G,CET],[CE,0]]
	double **K;
	K = (double**) malloc ((n+p)*sizeof(double*));
	// allocate and construct [G,CET]
	for (row=0; row<n; row++) {
		K[row] = (double*) malloc((n+p)*sizeof(double));
		for (col=0; col<n; col++)
			K[row][col] = G[row][col];
		for (col=n; col<n+p; col++)
			K[row][col] = CE[col-n][row];
	}
	// allocate and construct [CE,0]
	for (row=n; row<n+p; row++) {
		K[row] = (double*) malloc((n+p)*sizeof(double));
		for (col=0; col<n; col++)
			K[row][col] = CE[row-n][col];
		for (col=n; col<n+p; col++)
			K[row][col] = 0.0;
	}

	// allocate and construct k0 = [[g0],[ce0]]
	// lcqp: min 0.5 * x G x + g0 x
	// linear_solver: Gx = g0
	double *k0; // p+n
	k0 = (double*) malloc ((n+p)*sizeof(double));
	for (row=0; row<n; row++)
		k0[row] = -g0[row]; // minus explained above
	// lcqp: CE^T x + ce0 = 0
	// linear_solver: Ax = c
	for (row=n; row<n+p; row++)
		k0[row] = -ce0[row-n]; // minus explained above

	// allocate and construct xl = [[x],[lambda]]
	double *xl; // p+n
	xl = (double*) malloc ((n+p)*sizeof(double));

	// printf("calling linear solver; n+p = %d, K =\n",n+p);
	// print_sq_mat(n+p,K);
	// printf("\nk0 =\n");
	// print_vec(n+p,k0);
	// printf("\n");

	int steps = conjugate_solver ( n+p, K, k0, xl );
	// extract original variables x, toss lagrange variables
	for (row=0; row<n; row++)
		x[row] = xl[row];
	for (row=n; row<n+p; row++)
		lagrange[row-n] = -xl[row];

	// print results
	// printf("dut lcqp: converged in %d steps.\n", steps);

	// score = min 0.5 * x G x + g0 x
	double score = vec_mat_vec_mult(n,x,G)/2 + vec_vec_mult(n,g0,x);
	// printf("dut lcqp: score = %f\n", score);

	free (xl);
	free (k0);
	for (row=0; row<n+p; row++)
		free(K[row]);
	free (K);
	return score;
}