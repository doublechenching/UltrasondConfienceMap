/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSE_SOLVER_INTERFACE_H__
#define SPARSE_SOLVER_INTERFACE_H__

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

/** \brief	Interface for sparse solvers that take Eigen matrices
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 */

class SparseSolverInterface
{
public:
	/// Solver random walks system LuX=b, matrix returned as vector in column-major order
	virtual std::vector<double> solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label) = 0;
private:
};

#endif