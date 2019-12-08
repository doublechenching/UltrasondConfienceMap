#include "SparseSolverEigenLLT.h"
#include <Eigen/SparseCholesky>

std::vector<double> SparseSolverEigenLLT::solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label)
{
	SimplicialLLT<SparseMatrix<double> > solver;

	// Convert to dense format for solver compatibility - this might become a memory issue
	VectorXd b_dense = b;

	// LLT Decomposition
	solver.compute(A);

	/*if(solver.info()!=true) {
	  std::cout << "Decomposition failed!\n";
	}*/

	// Solve system with decomposition
	VectorXd x_dense = solver.solve(b_dense);

	/*if(solver.info()!=true) {
		std::cout << "Solver failed!\n";
	}*/

	// solve for another right hand side - good feature for multi-label setup
	//x1 = solver.solve(b1);

	// Holder for saving to external format
	std::vector<double> xmat(numel);

	for (int i=0; i<x_dense.rows(); i++)
	{
		double val = x_dense(i);
		xmat[uidx[i]] = val;
	}

	for (int i=0; i<seeds->size(); i++)
	{
		if((*labels)[i] == active_label)
			xmat[(*seeds)[i]] = 1.0;
		else
			xmat[(*seeds)[i]] = 0.0;
	}

	return xmat;
}