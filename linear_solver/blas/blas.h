#include <stdio.h>
#include <math.h>

void transpose (
	int n,
	double **A,
	double **AT
);

double vec_mat_vec_mult (
	int n,
	double *x,
	double **A
);

void mat_vec_mult (
	int n,
	double **A,
	double *x,
	double *result
);

double vec_vec_mult (
	int n,
	double *a,
	double *b
);

void mat_mat_mult (
	int n,
	double **A,
	double **B,
	double **result
);

void inverse (
	int n,
	double ** inp,
	double ** out
);

void inverse_2x2 (
	double ** inp,
	double ** out
);

void inverse_3x3 (
	double ** inp,
	double ** out
);

double two_norm (
	int n,
	double *a,
	double *b
);

void print_vec (
	int n,
	double *b
);

void print_mat (
	int n,
	int m,
	double **A
);

void print_sq_mat (
	int n,
	double **A
);