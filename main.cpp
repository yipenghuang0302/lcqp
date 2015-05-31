#include <stdio.h>
#include "lcqp.h"
// #include "quadprog/src/QuadProg++.hh"

int main () {

	srand (time(0));

	// Generate random Ax = b
	int n = 8; // size of A matrix
	int p = 8; // number of equality constraints
	int m = 8; // number of inequality constraints
	int row, col;
	double **A, **AT, **ATA, *b, **E, *e, **I, *i, *dut_x; // the A, aka G, matrix for quadprog has to be SPD

	A = (double**) malloc (n*sizeof(double*));
	AT = (double**) malloc (n*sizeof(double*));
	ATA = (double**) malloc (n*sizeof(double*));
	b = (double*) malloc (n*sizeof(double));
	E = (double**) malloc (p*sizeof(double*)); // p*n
	e = (double*) malloc (p*sizeof(double)); // p
	I = (double**) malloc (m*sizeof(double*)); // m*n
	i = (double*) malloc (m*sizeof(double)); // m
	dut_x = (double*) malloc (n*sizeof(double));

	// quadprog solver
	// QuadProgPP::Matrix<double> G, CE, CI;
	// QuadProgPP::Vector<double> g0, ce0, ci0, tb_x;
	// G.resize(n,n); CE.resize(n,p); CI.resize(n,m);
	// g0.resize(n); ce0.resize(p); ci0.resize(m);
	// tb_x.resize(n);

	// allocate and randomize
	for (row=0; row<n; row++) {
		A[row] = (double*) malloc(n*sizeof(double));
		AT[row] = (double*) malloc(n*sizeof(double));
		ATA[row] = (double*) malloc(n*sizeof(double));
		for (col=0; col<n; col++)
			A[row][col] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
		b[row] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
	}

	// allocate and randomize equality constraints
	for (row=0; row<p; row++) {
		E[row] = (double*) malloc(n*sizeof(double));
		for (col=0; col<n; col++) {
			E[row][col] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
			// CE[col][row] = E[row][col]; // the testbench expects transpose
		}
		e[row] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
		// ce0[row] = e[row];
	}

	// allocate and randomize inequality constraints
	for (row=0; row<m; row++) {
		I[row] = (double*) malloc(n*sizeof(double));
		for (col=0; col<n; col++) {
			I[row][col] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
			// CI[col][row] = I[row][col]; // the testbench expects transpose
		}
		i[row] = (double)(rand()-RAND_MAX/2) / (double)RAND_MAX;
		// ci0[row] = i[row];
	}

	transpose (n, A, AT);
	mat_mat_mult (n, AT, A, ATA);
	// copy to tb quadprog variables
	for (row=0; row<n; row++) {
		for (col=0; col<n; col++) {
			// G[row][col] = ATA[row][col];
			;
		}
		// g0[row] = b[row];
	}

	printf("problem: A =\n"); print_sq_mat(n,ATA); printf("\n");
	printf("problem: b = "); print_vec(n,b); printf("\n");
	printf("problem: E =\n"); print_mat(p,n,E); printf("\n");
	printf("problem: e = "); print_vec(p,e); printf("\n");
	printf("problem: I =\n"); print_mat(m,n,I); printf("\n");
	printf("problem: i = "); print_vec(m,i); printf("\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	/***********************
	** tb quadprog solver **
	************************/
	// std::cout << "tb quadprog: score = " << solve_quadprog(G, g0, CE, ce0, CI, ci0, tb_x) << std::endl;
	// std::cout << "tb quadprog: x = " << tb_x << std::endl;

	/* FOR DOUBLE CHECKING COST since in the solve_quadprog routine the matrix G is modified */
	// for (row=0; row<n; row++)
	// 	for (col=0; col<n; col++)
	// 		G[row][col] = ATA[row][col];

	// std::cout << "tb quadprog: score = ";
	// double sum = 0.0;
	// for (row=0; row<n; row++)
	// 	for (col=0; col<n; col++)
	// 		sum += tb_x[row] * G[row][col] * tb_x[col];
	// sum *= 0.5;	

	// for (row=0; row<n; row++)
	// 	sum += g0[row] * tb_x[row];
	// std::cout << sum << std::endl;
	// printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	/**********************
	*** dut lcqp solver ***
	***********************/
	lcqp (n, ATA, b, p, E, e, m, I, i, dut_x);

	// validate result
	mat_vec_mult(n, ATA, dut_x, b); // overwrites b
	printf("dut lcqp: x = "); print_vec(n,dut_x); printf("\n");
	printf("dut lcqp: b = "); print_vec(n,b); printf("\n");

	/************
	** cleanup **
	************/
	for (row=0; row<n; row++) {
		free(ATA[row]);
		free(AT[row]);
		free(A[row]);
	}
	for (row=0; row<p; row++)
		free(E[row]);
	for (row=0; row<m; row++)
		free(I[row]);
	free (dut_x);
	free (i);
	free (I);
	free (e);
	free (E);
	free (b);
	free (ATA);
	free (AT);
	free (A);
	return 0;

}