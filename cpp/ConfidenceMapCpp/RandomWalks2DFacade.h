/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef RANDOM_WALKS_2D_FACADE_H__
#define RANDOM_WALKS_2D_FACADE_H__

/** \brief	Compute random walks probabilities for 2D image
 *	\author	Athanasios Karamalis
 *  \date	10.06.2012
 *
 *	Implementing the paper of
 *	Grady, L.: Random walks for image segmentation, IEEE Transactions on Pattern Analysis and Machine Intelligence, 28, 11, 1768--1783, 2006
 */
#include <string>
#include <vector>
class RandomWalks2D;

class RandomWalks2DFacade
{
public:
	RandomWalks2DFacade();
	virtual ~RandomWalks2DFacade();

	/// Compute random walks probabilities for image
	std::vector<double> computeProbabilities(int active_label, double beta= 90 );

	/// Set external image 2D as vector with column-major ordering
	void setImage(std::vector<double> &image, int rows, int cols);

	/// Set Seeds, labels, active label and #labels
	void setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels);

	/// Set desired solver (default Eigen LLT)
	void setSolver(std::string solver, int iterations = 2000, double tolerance = 1.0e-7);
private:
	RandomWalks2D * walks2D; ///< Implementation of random walks
};

#endif