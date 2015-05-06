#include "linear_solver.h"

int conjugate_solver (
	int n,
	double ** A,
	double * b,
	double * x
) {

	int row, col;
	double **AT, *ATb, **ATA;
	double *residual, *dee;

	AT = (double**) malloc (n*sizeof(double*));
	ATb = (double*) malloc (n*sizeof(double));
	ATA = (double**) malloc (n*sizeof(double*));
	residual = (double*) malloc (n*sizeof(double));
	dee = (double*) malloc (n*sizeof(double));
	
	// allocate
	for (row=0; row<n; row++) {
		AT[row] = (double*) malloc(n*sizeof(double));
		ATA[row] = (double*) malloc(n*sizeof(double));
	}

	// aux variables
	transpose (n, A, AT);
	mat_mat_mult (n, AT, A, ATA);
	mat_vec_mult (n, AT, b, ATb);

	// x is the vector that minimizes f(x) = 0.5(x^T)ATAx - ((ATb)^T)x
	// initialize
	int step = 0;
	for (row=0; row<n; row++) {
		residual[row] = ATb[row];
		for (col=0; col<n; col++)
			residual[row] -= ATA[row][col]*x[col]; // in HCDC, x elements are non constant all but row number of terms are constant
		dee[row] = residual[row];
	}

	// delta = rTr;
	double delta_new = vec_vec_mult (n,residual,residual);
	double delta_0 = delta_new;

	do {
		double alpha = delta_new / vec_mat_vec_mult (n, dee, ATA);// + DBL_EPSILON;
		for (row=0; row<n; row++)
			x[row] += alpha * dee[row];
		for (row=0; row<n; row++) {
			residual[row] = ATb[row];
			for (col=0; col<n; col++)
				residual[row] -= ATA[row][col]*x[col]; // in HCDC, x elements are non constant all but row number of terms are constant
		}
		double delta_old = delta_new;
		delta_new = vec_vec_mult (n,residual,residual);
		double beta = delta_new / delta_old;
		for (row=0; row<n; row++) {
			dee[row] = residual[row] + beta*dee[row];
		}
		step++;
		// printf("step = %d\n",step++);
	} while (
		delta_new>DBL_EPSILON*delta_0
		&& step<(1<<16)
	);

	for (row=0; row<n; row++) {
		free(AT[row]);
		free(ATA[row]);
	}
	free (AT);
	free (ATA);
	free (ATb);
	free (residual);
	free (dee);

	return step;

}