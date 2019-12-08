/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CONFIDENCE_MAPS_2D_FACADE_H__
#define CONFIDENCE_MAPS_2D_FACADE_H__
#include <string>
#include <vector>
/** \brief	Facade for computing confidence map for 2D image
 *	\author	Athanasios Karamalis
 *  \date	11.06.2012
 *
 *	Implementing ultrasound confidence estimation algorithm as described in
 *	A. Karamalis, W. Wein, T. Klein, N. Navab
 *	Ultrasound Confidence Maps using Random Walks
 *	Medical Image Analysis, 16, 6, 1101 - 1112, 2012, DOI: http://dx.doi.org/10.1016/j.media.2012.07.005 (bib)
 *
 *	Use the facade to automatically generate the boundary conditions for the confidence estimation problem.
 */

class ConfidenceMaps2D;

class ConfidenceMaps2DFacade
{
public:
	ConfidenceMaps2DFacade();
	virtual ~ConfidenceMaps2DFacade();

	/// Compute map for image
	std::vector<double> computeMap(double beta = 100, double gamma = 0.06);

	/// Set external image 2D as vector with column-major ordering
	void setImage(std::vector<double> &image, int rows, int cols, double alpha=2.0, bool normalizeValues=false );

	/// Set desired solver (default Eigen LLT)
	void setSolver(std::string solver, int iterations = 2000, double tolerance = 1.0e-7);
private:
	ConfidenceMaps2D * maps2D; ///< Implementation of confidence estimation
	std::vector<int> seeds; ///< Seeds/boundary conditions
	std::vector<int> labels; ///< Labels for seeds
};

#endif