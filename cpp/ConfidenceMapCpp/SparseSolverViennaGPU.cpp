#include "SparseSolverViennaGPU.h"

#define VIENNACL_HAVE_EIGEN
#include <viennacl/scalar.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/cg.hpp>

SparseSolverViennaGPU::SparseSolverViennaGPU(int iterations, double tolerance)
{
	this->iterations = iterations;
	this->tolerance = tolerance;
}

std::vector<double> SparseSolverViennaGPU::solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label)
{
	SparseMatrix<double> A_sparse_matrix(A);
	VectorXd b_dense = b;

	viennacl::linalg::cg_tag custom_tag(tolerance, iterations);

	viennacl::compressed_matrix<double>     viennacl_A(A_sparse_matrix.rows(), A_sparse_matrix.cols());
	viennacl::vector<double>                viennacl_b(b_dense.rows());
	viennacl::vector<double>                viennacl_x(b_dense.rows());

	viennacl::copy(A_sparse_matrix, viennacl_A);
	viennacl::copy(b_dense, viennacl_b);

	viennacl_x = viennacl::linalg::solve(viennacl_A, viennacl_b, custom_tag); // GPU Solver

	VectorXd x_dense(b.rows());
	viennacl::copy(viennacl_x,  x_dense);

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