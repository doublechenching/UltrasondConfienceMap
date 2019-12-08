/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CONFIDENCE_MAPS_2D_H__
#define CONFIDENCE_MAPS_2D_H__

#include "RandomWalksCore.h"

/** \brief	Confidence Map for 2D image
 *	\author	Athanasios Karamalis
 *  \date	11.06.2012
 *
 *	Implementing ultrasound confidence estimation algorithm as described in
 *	A. Karamalis, W. Wein, T. Klein, N. Navab
 *	Ultrasound Confidence Maps using Random Walks
 *	Medical Image Analysis, 16, 6, 1101 - 1112, 2012, DOI: http://dx.doi.org/10.1016/j.media.2012.07.005 (bib)
 */

class ConfidenceMaps2D : public RandomWalksCore
{
public:
	ConfidenceMaps2D();
	virtual ~ConfidenceMaps2D();

	/// Set data matrix. We assume normalized values. If this is not the case use the flag in the function call.
	void setMatrix2D(const std::vector<double> * matrix, int rows, int cols, bool normalizeValues=false);

	/// Implement the random walks 2D edge assembly (L is then finalized in RandomWalksCore)
	virtual void assembleLaplacianEdges();

	/// Set the alpha parameter - penalizing vertical walks with distance
	void setAlpha(double alpha){this->alpha = alpha;};

	/// Set the beta parameter (default random walks)
	void setBeta(double beta){this->beta = beta;};

	/// Set the gamma parameter - penalizing horizontal walks
	void setGamma(double gamma){this->gamma = gamma;};
private:
	double alpha; ///< Vertical walks penalty
	double beta; ///< Default walks parameter
	double gamma; ///< Horizontal walks penalty
};

#endif