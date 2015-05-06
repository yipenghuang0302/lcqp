#include <stdio.h>
#include <math.h>
#include "linear_solver/linear_solver.h"

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
);

double eq_con_lcqp (
	int n,
	double ** G,
	double * g0,
	int p,
	double ** CE,
	double * ce0,
	double * x,
	double * lagrange
);

int working_set (
	int n,
	int p,
	double ** CE,
	double * ce0,
	int m,
	double ** CI,
	double * ci0,
	unsigned int working_flags,
	double ** CECI,
	double * ce0ci0
);