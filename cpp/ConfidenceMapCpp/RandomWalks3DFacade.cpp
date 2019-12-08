#include "RandomWalks3DFacade.h"
#include "RandomWalks3D.h"

RandomWalks3DFacade::RandomWalks3DFacade()
{
	walks3D = new RandomWalks3D;
}

RandomWalks3DFacade::~RandomWalks3DFacade()
{
	if(walks3D!=0)
		delete walks3D;
	walks3D = 0;
}

void RandomWalks3DFacade::setSolver(std::string solver, int iterations, double tolerance)
{
	walks3D->setSolver(solver,iterations,tolerance);
}

std::vector<double> RandomWalks3DFacade::computeProbabilities(int active_label, double beta)
{
	walks3D->setActiveLabel(active_label);
	walks3D->setBeta(beta);

	std::vector<double> probabilities;

	probabilities = walks3D->solve();

	return probabilities;
}

void RandomWalks3DFacade::setImage(std::vector<double> &image, int rows, int cols, int stacks)
{
	walks3D->setMatrix3D(&image,rows,cols,stacks);
}

void RandomWalks3DFacade::setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels)
{
	walks3D->setLabeling(seeds, labels, active_label, num_labels);
}