#include "lcqp.h"

/*
	The problem is in the form:
	min 0.5 * xT G x + g0 x
	(1/2) xT P x + qT x + r
	s.t.
	CE x + ce0 = 0
	CI x + ci0 >= 0

	The matrix and vectors dimensions are as follows:
	G: n*n
	g0: n
	CE: p*n
	ce0: p
	CI: m*n
	ci0: m
	x: n
*/

double lcqp (
	int n,
	double ** G,
	double * g0,
	int p,
	double ** CE,
	double * ce0,
	int m,
	double ** CI,
	double * ci0,
	double * x
) {

	int row;

	// create zero vector of size n
	double *zeron = (double*) malloc (n*sizeof(double));
	for (row=0; row<n; row++)
		zeron[row] = 0.0;
	// create zero vector of size m+p
	double *zeromp = (double*) malloc ((m+p)*sizeof(double));
	for (row=0; row<m+p; row++)
		zeromp[row] = 0.0;

	// allocate CECI, the working set of constraints
	// all equality plus some random inequalities treated as equality constraints
	double **CECI = (double**) malloc ((m+p)*sizeof(double*)); // (m+p)*n
	for (row=0; row<m+p; row++)
		CECI[row] = (double*) malloc (n*sizeof(double));
	double *ce0ci0 = (double*) malloc ((m+p)*sizeof(double)); // m+p

	// Compute a feasible starting point x0;
	// Set W0 to be a subset of the active constraints at x0;
	int working_flags = 0;
	for (row=0; row<n-p; row++)
		working_flags += 1<<row;
	working_set ( n, p, CE, ce0, m, CI, ci0, working_flags, CECI, ce0ci0 );

	// initialize with optimum of equality constrained optimization
	double *lagrange = (double*) malloc ((m+p)*sizeof(double));
	double score = eq_con_lcqp ( n, G, g0, m+p, CECI, ce0ci0, x, lagrange );

	double *g_k = (double*) malloc (n*sizeof(double));
	double *p_k = (double*) malloc (n*sizeof(double));
	bool done = false;
	do {

		/**find step size**/
		// g_k = Gx_k + c
		mat_vec_mult ( n, G, x, g_k );
		for (row=0; row<n; row++)
			g_k[row] += g0[row];
		// optimize
		score = eq_con_lcqp ( n, G, g_k, m+p, CECI, zeromp, p_k, lagrange );

		/**if p_k=0.0, no steps to take**/
		if (two_norm(n,p_k,zeron)<FLT_EPSILON) {

			// see if there are inequality constraints in working set to remove
			// determine the smallest (most negative) lagrange multiplier
			int least_lagrange_index = 0;
			double least_lagrange_value = 0.0;
			for (row=0; row<m; row++)
				if ((working_flags & 1<<row) != 0) // if inequality constraint active
					if (lagrange[row+p] < least_lagrange_value) {
						least_lagrange_index = row;
						least_lagrange_value = lagrange[row+p];
					}

			// remove the constraint, reiterate and find step size again
			if (least_lagrange_value < 0.0) {
				// remove the most negative lagrange value from constraints
				working_flags -= 1<<least_lagrange_index;
				working_set ( n, p, CE, ce0, m, CI, ci0, working_flags, CECI, ce0ci0 );
			} else {
				// all lagranges non-negative, so optimum is with current working set
				done = true;
			}

		/**there is non zero step to take**/
		} else {

			// find alpha, the fraction of the step size allowed to take
			// consider inequality constraints not in working set
			int blocking_index = 0;
			double least_alpha = 1.0;
			for (row=0; row<m; row++)
				if ((working_flags & 1<<row)==0 && vec_vec_mult(n,CI[row],p_k)<0.0) {
					double step_fraction = (double) (-ci0[row]-vec_vec_mult(n,CI[row],x)) / (double) vec_vec_mult(n,CI[row],p_k);
					if (step_fraction < least_alpha) {
						blocking_index = row;
						least_alpha = step_fraction;
					}
				}

			// now alpha is most conservative step size fraction
			for (row=0; row<n; row++) {
				// printf("least_alpha = %f\n", least_alpha);
				x[row] += least_alpha * p_k[row];
			}

			// add the constraint, reiterate and find step size again
			if (least_alpha < 1.0) {
				// add the most blocking inequality to constraints
				working_flags += 1<<blocking_index;
				working_set ( n, p, CE, ce0, m, CI, ci0, working_flags, CECI, ce0ci0 );
			}
		}

	} while (!done);

	free (lagrange);
	free (p_k);
	free (g_k);

	free (ce0ci0);
	for (row=0; row<m+p; row++)
		free (CECI[row]);
	free (CECI);

	free (zeromp);
	free (zeron);

	return score;
}