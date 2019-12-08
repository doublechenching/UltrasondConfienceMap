/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSE_SOLVER_BIGSTAB_H__
#define SPARSE_SOLVER_BIGSTAB_H__

#include "SparseSolverInterface.h"

/** \brief	Bi-conjugate gradient stabilized solver for sparse square problems. Eigen solver for random walks system
 *	\author	Athanasios Karamalis
 *  \date	17.07.2012
 */

class SparseSolverEigenBiCGSTAB : public SparseSolverInterface
{
public:
	SparseSolverEigenBiCGSTAB(int iterations, double tolerance);
	/// Solver with BiGSTAB
	virtual std::vector<double> solve_Ax_b(SparseMatrix<double> A, SparseVector<double> b, int numel, std::vector<int> & uidx, const std::vector<int> * labels, const std::vector<int> * seeds, int active_label);
protected:

private:
	int iterations; ///< CG iterations
	double tolerance; ///< CG tolerance
};

#endif