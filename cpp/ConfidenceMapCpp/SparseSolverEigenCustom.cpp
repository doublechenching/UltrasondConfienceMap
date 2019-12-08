#include "SparseSolverEigenCustom.h"

SparseSolverEigenCustom::SparseSolverEigenCustom(int iterations, double tolerance)
{
	this->iterations = iterations;
	this->tolerance = tolerance;
}

std::vector<double> SparseSolverEigenCustom::solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label)
{
	// Custom solver
	SparseVector<double> x = this->cgSolveSparse(A,b,iterations,tolerance);

	std::vector<double> xmat(numel);

	// All entries are non-zero except at seeds
	for (SparseVector<double>::InnerIterator it(x); it; ++it)
	{
			xmat[uidx[it.index()]] = it.value();
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

SparseVector<double> SparseSolverEigenCustom::cgSolveSparse(const SparseMatrix<double> & A,const SparseVector<double> & b,int iter, double residual)
{
	SparseVector<double> r(b.rows());
	SparseVector<double> p(b.rows());
	SparseVector<double> Ap(b.rows());
	SparseVector<double> x(b.rows());

	r = b - A *x;
	p = r;

	double rTr,pTAp,alpha,beta,rTrnew,rnorm;
	SparseVector<double> vtemp;
	bool isConverged = false;
	for(int k=0;k<iter;k++)
	{
		Ap = A*p;
		vtemp = r.transpose()*r;
		rTr = vtemp.coeff(0);

		vtemp = p.transpose()*Ap;
		pTAp = vtemp.coeff(0);
		alpha = rTr/pTAp;

		x = x + (alpha * p);
		r = r - (alpha * Ap);
		rnorm = r.norm();
		if(rnorm<residual)
		{
			isConverged = true;
			break;
		}

		vtemp = r.transpose()*r;
		rTrnew = vtemp.coeff(0);

		beta = rTrnew / rTr;
		p = r + (beta * p);
	}

	return x;
}