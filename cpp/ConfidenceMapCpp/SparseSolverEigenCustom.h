/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSE_SOLVER_EIGEN_CUSTOM_H__
#define SPARSE_SOLVER_EIGEN_CUSTOM_H__

#include "SparseSolverInterface.h"

/** \brief	Conjugate Gradient Custom Eigen solver for random walks system
 *	\author	Athanasios Karamalis
 *  \date	11.06.2012
 */

class SparseSolverEigenCustom : public SparseSolverInterface
{
public:
	SparseSolverEigenCustom(int iterations, double tolerance);
	/// Solver with CG
	virtual std::vector<double> solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label);
private:
	SparseVector<double> cgSolveSparse(const SparseMatrix<double> & A,const SparseVector<double> & b,int iter, double residual);
	int iterations; ///< CG iterations
	double tolerance; ///< CG tolerance
};

#endif