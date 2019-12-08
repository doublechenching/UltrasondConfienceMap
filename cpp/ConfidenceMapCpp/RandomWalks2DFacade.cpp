#include "RandomWalks2DFacade.h"
#include "RandomWalks2D.h"

RandomWalks2DFacade::RandomWalks2DFacade()
{
	walks2D = new RandomWalks2D;
}

RandomWalks2DFacade::~RandomWalks2DFacade()
{
	if(walks2D!=0)
		delete walks2D;
	walks2D = 0;
}

void RandomWalks2DFacade::setSolver(std::string solver, int iterations, double tolerance)
{
	walks2D->setSolver(solver,iterations,tolerance);
}

std::vector<double> RandomWalks2DFacade::computeProbabilities(int active_label, double beta)
{
	walks2D->setActiveLabel(active_label);
	walks2D->setBeta(beta);

	std::vector<double> probabilities;

	probabilities = walks2D->solve();

	return probabilities;
}

void RandomWalks2DFacade::setImage(std::vector<double> &image, int rows, int cols)
{
	walks2D->setMatrix2D(&image,rows,cols);
}

void RandomWalks2DFacade::setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels)
{
	walks2D->setLabeling(seeds,labels,active_label,num_labels);
}