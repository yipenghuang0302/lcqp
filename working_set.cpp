#include "lcqp.h"

int working_set (
	int n,
	int p,
	double ** CE, // p*n
	double * ce0, // p
	int m,
	double ** CI, // m*n
	double * ci0, // m
	unsigned int working_flags,
	double ** CECI, // (p+m)*n
	double * ce0ci0 // p+m
) {

	int row, col, working_count=0;

	// copy the equality constraint set
	for (row=0; row<p; row++)
		for (col=0; col<n; col++) {
			CECI[row][col] = CE[row][col];
			ce0ci0[row] = ce0[row];
		}

	// copy the inequality constraint set
	for (row=0; row<m; row++)
		// working
		if ((working_flags & 1<<row) != 0) {
			working_count++;
			for (col=0; col<n; col++)
				CECI[row+p][col] = CI[row][col];
			ce0ci0[row+p] = ci0[row];
		// not working
		} else {
			for (col=0; col<n; col++)
				CECI[row+p][col] = 0.0;
			ce0ci0[row+p] = 0.0;
		}

	printf ("working_flags = %d; set size = %d\n", working_flags, p+working_count);
	return p+working_count;
}