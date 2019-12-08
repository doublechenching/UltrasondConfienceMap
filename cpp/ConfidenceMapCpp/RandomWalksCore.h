/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef RANDOM_WALKS_CORE_H__
#define RANDOM_WALKS_CORE_H__

#ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif

#include <Eigen/Sparse>
#include <vector>
#include <string>
#include "SparseSolverInterface.h"

/** \brief	Core functionality for creating the Graph-Laplacian for the random walks system
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 *
 *	Implementation based on the work:
 *	Grady, L.: Random walks for image segmentation, IEEE Transactions on Pattern Analysis and
 *	Machine Intelligence,28,11,1768-1783,2006
 *
 *	Performs most of the random walks algorithmic operations, except of assembling the Laplacian matrix.
 *	The Laplacian is problem specific and dimension specific.
 *
 *	Different solvers can be used for Ax=b including Conjugate Gradient for Eigen, the CG solver from ViennaCL
 *	that is available for both CPU and GPU
 */
using namespace Eigen;

class SparseSolverFactory;
class RandomWalksCore
{
public:
	RandomWalksCore();
	virtual ~RandomWalksCore();
	/// Choose between different solvers
	/**
		Eigen-LLT
		Eigen-CG
		Eigen-CG-Custom
		Vienna-CG-CPU
		Vienna-CG-GPU
	*/
	void setSolver(std::string solver, int iterations = 2000, double tolerance = 1.0e-7);

	/// Solve the random walks problem, i.e, Lu x = -B^T M
	std::vector<double> solve();

	/// Set labeling including seeds
	/** The matrix is given in each concrete implementation as it can be 2D and 3D */
	void setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels);

	/// Reset active label (solving for multiple labels)
	void setActiveLabel(int active_label){this->active_label = active_label;};
protected:
	/// Virtual function for assembling the graph Laplacian edges.
	/** Here you can define the graph connectivity and the edge weighting. */
	virtual void assembleLaplacianEdges() = 0;

	/// Compute and fill degree in Laplacian matrix
	/** After inserting the edge connections and weights in the Laplacian simply compute the degree and fill the diagonal with it */
	void assembleDegree();

	/// Assemble block matrix Lu and the rhs for the solver b.
	/*
	Sorted Laplace matrix for random walks/circuit problem

	| Lm(nxn)	B(nxq)  |
	| B^T(qxn)	Lu(qxq) |

	Lu/Lm : relation of un-/marked nodes
	n : # marked nodes
	q : # unmarked nodes

	The final solution to the random walks problem is given as
	Lu X = -B^T M
	where X(qxl) are the potentials at the unmarked nodes with q is #unmarked nodes and l is #labels,
	M(nxl) are the marked unit potentials. For l=2 M(n,1) = 1 for n marked node with label 1.

	Instead of re-ordering the Laplace matrix (L) to get the blocks we can get them individually. For example, to get B^T we assemble it by taking
	the unmarked rows and the marked columns of L. For Lm we take the marked rows and marked cols from L and for Lu the unmarked rows and unmarked cols from L.

	Accessing entries in sparse matrices in Eigen using for example coeff costs log(rho*outer_size)!

	Therefore, we assemble the matrices through matrix operations with permutation matrix.
	*/
	void assemble_Lu_b();

	/// Generates which indices are unmarked
	void generateUniqueIndices();

	SparseMatrix<double> L; ///< Laplacian matrix

	const std::vector<double> * matrix; ///< Pointer to matrix containing the data - column major assumed

	int rows; ///< # rows

	int cols; ///< # cols

	int numel; ///< # numel = rows * cols

private:

	SparseMatrix<double> Lu; ///< Block Matrix Lu from L

	SparseVector<double> b; ///< rhs solution

	SparseSolverInterface * sparseSolver; ///< Sparse solver to be used for Ax=b

	/// Pointer to seeds. Seeds assumed to be unique and positive.
	/** They are indices to the matrix data/image entries assuming column-major order. */
	const std::vector<int> * seeds;

	/// Labels, i.e., for each seed element its corresponding label is given. Thus seeds.size = labels.size
	const std::vector<int> * labels;

	/// Currently active label for which we get the random walks solution.
	/** Future implementation should allow an LU decomposition for the Lu matrix so that the
	* solution can be used for any label with simple matrix-vector multiplication instead of solving
	* the entire system of equations again.
	*/
	int active_label;

	int num_labels; ///< Number of unique labels

	std::vector<int> uidx; ///< Matrix indices for unmarked nodes

	/// Is the Laplacian and Lu available
	/** Once the Laplacian for given seeds is assembled we can use it for solving for different labels */
	bool isLaplaceAvailable;

	SparseSolverFactory * solverFactory; ///< Factory for different solvers
};

#endif