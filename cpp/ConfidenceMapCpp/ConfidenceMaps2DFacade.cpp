#include "ConfidenceMaps2DFacade.h"
#include "ConfidenceMaps2D.h"

ConfidenceMaps2DFacade::ConfidenceMaps2DFacade()
{
	maps2D = new ConfidenceMaps2D();
}

ConfidenceMaps2DFacade::~ConfidenceMaps2DFacade()
{
	if(maps2D!=0)
		delete maps2D;
	maps2D = 0;
}

std::vector<double> ConfidenceMaps2DFacade::computeMap(double beta, double gamma)
{
	maps2D->setBeta(beta);
	maps2D->setGamma(gamma);

	std::vector<double> map;
	map = maps2D->solve();
	return map;
}

void ConfidenceMaps2DFacade::setImage(std::vector<double> &image, int rows, int cols, double alpha, bool normalizeValues)
{
	maps2D->setAlpha(alpha);
	maps2D->setMatrix2D(&image,rows,cols,normalizeValues);
	this->seeds.clear();
	this->labels.clear();

	// Automatic seed generation for confidence maps

	seeds.reserve(cols*2);
	labels.reserve(cols*2);

	// Upper row seeds
	for(int i=0;i<cols;i++)
	{
		seeds.push_back(i*rows);
	}

	// Lower row seeds
	for(int i=0;i<cols;i++)
	{
		seeds.push_back((i*rows)+(rows-1));
	}

	// UpperLabels
	for(int i=0;i<cols;i++)
	{
		labels.push_back(0);
	}

	// Lover labels
	for(int i=0;i<cols;i++)
	{
		labels.push_back(1);
	}

	maps2D->setLabeling(&seeds,&labels,0,2);
}

void ConfidenceMaps2DFacade::setSolver(std::string solver, int iterations, double tolerance)
{
	maps2D->setSolver(solver,iterations,tolerance);
}