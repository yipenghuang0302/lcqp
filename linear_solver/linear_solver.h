#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include "blas/blas.h"

int conjugate_solver (
	int n,
	double ** A,
	double * b,
	double * x_next
);

int steepest_solver (
	int n,
	double ** A,
	double * b,
	double * x_next
);