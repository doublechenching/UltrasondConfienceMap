/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef RANDOM_WALKS_3D_H__
#define RANDOM_WALKS_3D_H__

#include <Eigen/Sparse>
#include "RandomWalksCore.h"

/** \brief	Random walks for 3D data matrix.
 *	\author	Athanasios Karamalis
 *  \date	11.06.2012
 *
 *	Assemple 6-connected Laplacian. Let base-class handles the rest.
 *
 *	Implementing the paper of
 *	Grady, L.: Random walks for image segmentation, IEEE Transactions on Pattern Analysis and Machine Intelligence, 28, 11, 1768--1783, 2006
 */

using namespace Eigen;

class RandomWalks3D : public RandomWalksCore
{
public:
	/// Set data matrix
	void setMatrix3D(const std::vector<double> * matrix, int rows, int cols, int stacks);

	/// Implement the random walks 2D edge assembly (L is then finalized in RandomWalksCore)
	virtual void assembleLaplacianEdges();

	/// Set the beta parameter
	void setBeta(double beta){this->beta = beta;};
private:
	double beta; ///< Default walks parameter

	int stacks; ///< Number of stacks/slices in 3D data matrix/volume

	int numel_stack; ///< Number of elements per stack
};

#endif