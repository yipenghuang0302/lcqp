#include "blas.h"

void transpose ( int n, double **A, double **AT ) {
	int row, col;
	for (row=0; row<n; row++)
		for (col=0; col<n; col++)
			AT[row][col] = A[col][row];
}

double vec_mat_vec_mult ( int n, double *x, double **A ) {
	double sum = 0.0;
	int row, col;
	for (row=0; row<n; row++)
		for (col=0; col<n; col++)
			sum += x[row] * A[row][col] * x[col];
	return sum;
}

void mat_vec_mult ( int n, double **A, double *x, double *result ) {
	int row, col;
	for (row=0; row<n; row++) {
		result[row] = 0.0;
		for (col=0; col<n; col++)
			result[row] += A[row][col] * x[col];
	}
}

double vec_vec_mult ( int n, double *a, double *b ) {
	double sum = 0.0;
	int row;
	for (row=0; row<n; row++)
		sum += a[row] * b[row];
	return sum;
}

void mat_mat_mult ( int n, double **A, double **B, double **result ) {
	int row, col, iter;
	for (row=0; row<n; row++)
		for (col=0; col<n; col++) {
			result[row][col] = 0.0;
			for (iter=0; iter<n; iter++)
				result[row][col] += A[row][iter] * B[iter][col];
		}
}

void inverse ( int n, double ** inp, double ** out ) {
	switch (n) {
		case 2: inverse_2x2 (inp, out); break;
		case 3: inverse_3x3 (inp, out); break;
		default: break;
	}
}

void inverse_2x2 ( double ** inp, double ** out ) {
	double ad = (inp[0][0])*(inp[1][1]);
	double bc = (inp[0][1])*(inp[1][0]);
	double det = ad-bc;
	double invdet = 1 / ( det );
	out[0][0] = (inp[1][1])*(invdet);
	out[0][1] = -(inp[0][1])*(invdet);
	out[1][0] = -(inp[1][0])*(invdet);
	out[1][1] = (inp[0][0])*(invdet);
}

void inverse_3x3 (double ** a, double ** out) {
	int i, j;
	double det = 0.0;
	for (i=0;i<3;i++) {
		det = det + (
			a[0][i] * (
				a[1][(i+1)%3] * a[2][(i+2)%3] -
				a[1][(i+2)%3] * a[2][(i+1)%3]
			)
		);
	}
	for (i=0;i<3;i++) {
		for (j=0;j<3;j++)
			out [j][i] = (
				a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3] - 
				a[(i+1)%3][(j+2)%3] * a[(i+2)%3][(j+1)%3]
			) / det;
	}
}

double two_norm ( int n, double *a, double *b ) {
	double sum = 0.0;
	int row;
	for ( row=0; row<n; row++ )
		sum += (a[row]-b[row]) * (a[row]-b[row]);
	return sqrt(sum);
}

void print_vec ( int n, double *b ) {
	if (n>0) printf("[%f", b[0]);
	int row;
	for (row=1; row<n; row++) printf(", %f", b[row]);
	printf("]");
}

void print_mat ( int n, int m, double **A ) {
	printf("[");
	if (n>0) print_vec(m, A[0]);
	int row;
	for (row=1; row<n; row++) {
		printf(",\n");
		print_vec(m, A[row]);
	}
	printf("]");
}

void print_sq_mat ( int n, double **A ) {
	printf("[");
	if (n>0) print_vec(n, A[0]);
	int row;
	for (row=1; row<n; row++) {
		printf(",\n");
		print_vec(n, A[row]);
	}
	printf("]");
}