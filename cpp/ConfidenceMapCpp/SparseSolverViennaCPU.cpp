#include "SparseSolverViennaCPU.h"

#define VIENNACL_HAVE_EIGEN
#include <viennacl/scalar.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>

SparseSolverViennaCPU::SparseSolverViennaCPU(int iterations, double tolerance)
{
	this->iterations = iterations;
	this->tolerance = tolerance;
}

std::vector<double> SparseSolverViennaCPU::solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label)
{
	SparseMatrix<double> A_sparse_matrix(A);
	VectorXd b_dense = b;

	viennacl::linalg::cg_tag custom_tag(tolerance, iterations);
	VectorXd x_dense = viennacl::linalg::solve(A_sparse_matrix, b_dense, custom_tag); // CPU Solver with dense

	std::vector<double> xmat(numel);

	for (int i=0; i<x_dense.rows(); i++)
	{
		xmat[uidx[i]] = x_dense(i);
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