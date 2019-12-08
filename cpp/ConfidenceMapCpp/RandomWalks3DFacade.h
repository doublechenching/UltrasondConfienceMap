/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef RANDOM_WALKS_3D_FACADE_H__
#define RANDOM_WALKS_3D_FACADE_H__

/** \brief	Compute random walks probabilities for 3D image
 *	\author	Athanasios Karamalis
 *  \date	11.06.2012
 */
#include <string>
#include <vector>
class RandomWalks3D;

class RandomWalks3DFacade
{
public:
	RandomWalks3DFacade();
	virtual ~RandomWalks3DFacade();

	/// Compute random walks probabilities for image
	std::vector<double> computeProbabilities(int active_label, double beta = 90);

	/// Set external image 3D as vector with column-major ordering
	void setImage(std::vector<double> &image, int rows, int cols, int stacks);

	/// Set Seeds, labels, active label and #labels
	void setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels);

	/// Set desired solver (default Eigen LLT)
	void setSolver(std::string solver, int iterations = 2000, double tolerance = 1.0e-7);
private:
	RandomWalks3D * walks3D; ///< Implementation of random walks
};

#endif