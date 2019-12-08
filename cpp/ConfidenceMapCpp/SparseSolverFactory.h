/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSE_SOLVER_FACTORY_H__
#define SPARSE_SOLVER_FACTORY_H__

#include "SparseSolverInterface.h"
#include <string>

/** \brief	Solver factory for instantiating the desired solver (Eigen, Vienna CPU/GPU).
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 */

class SparseSolverFactory
{
public:
	// Choose between different solvers
	/**
		Eigen-LLT
		Eigen-CG
		Eigen-BiCGSTAB
		Eigen-CG-Custom
		Vienna-CG-CPU
		Vienna-CG-GPU
	*/
	SparseSolverInterface * createSolver(std::string type="Eigen-LLT", int iterations = 3000, double tolerance = 1.0e-7);

	/// Set parameters for conjugate gradient solvers
	void setConjugateGradientParams(int iterations, double tolerance);

	/// Get iterative solver number of iterations
	int getIterations(){return iterations;};

	/// Get iterative solver tolerance
	double getTolerance(){return tolerance;};
private:

	int iterations; ///< Iterations for CG solvers

	double tolerance; ///< Tolerance for CG solvers
};

#endif