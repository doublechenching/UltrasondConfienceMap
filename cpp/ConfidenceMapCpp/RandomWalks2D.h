/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef RANDOM_WALKS_2D_H__
#define RANDOM_WALKS_2D_H__

#include <Eigen/Sparse>
#include "RandomWalksCore.h"

/** \brief	Random walks for 2D data matrix.
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 *
 *	Note: Assembling of 2D graph-laplacian has been trible-checked and the Lu and b
 *	matrices are exactly the same as in Matlab up to 1.0e-14 precision.
 *
 *	Assemple 4-connected Laplacian. Let base-class handle the rest.
 *
 *	Implementing the paper of
 *	Grady, L.: Random walks for image segmentation, IEEE Transactions on Pattern Analysis and Machine Intelligence, 28, 11, 1768--1783, 2006
 */

using namespace Eigen;

class RandomWalks2D : public RandomWalksCore
{
public:
	/// Set data matrix in column-major ordering.
	void setMatrix2D(const std::vector<double> * matrix, int rows, int cols);

	/// Implement the random walks 2D edge assembly (L is then finalized in RandomWalksCore)
	virtual void assembleLaplacianEdges();

	/// Set the beta parameter
	void setBeta(double beta){this->beta = beta;};
private:
	double beta;	///< Default beta parameter
};

#endif