/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSE_SOLVER_EIGEN_LLT_H__
#define SPARSE_SOLVER_EIGEN_LLT_H__

#include "SparseSolverInterface.h"

/** \brief	Direct Eigen solver using LLT  for random walks system
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 */

class SparseSolverEigenLLT : public SparseSolverInterface
{
public:
	/// Solver with direct LLT (Cholesky Decomposition)
	// TODO Save decomposition for next solution
	virtual std::vector<double> solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label);
};
#endif