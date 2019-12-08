#include "RandomWalks3D.h"

void RandomWalks3D::setMatrix3D(const std::vector<double> * matrix, int rows, int cols, int stacks)
{
	this->matrix = matrix;
	this->rows = rows;
	this->cols = cols;
	this->stacks = stacks;
	this->numel = rows*cols*stacks;
	this->numel_stack = rows*cols;
}

void RandomWalks3D::assembleLaplacianEdges()
{
	//const double epsilon = std::numeric_limits<double>::epsilon(); // Numerical limit to avoid division by zero and add to the weights to avoid zero weights in Laplacian
	// Constant from Matlab
	const double epsilon = 1e-5; // Numerical limit to avoid division by zero and add to the weights to avoid zero weights in Laplacian

	this->L.resize(numel,numel);
	double min_weight = 0;
	double max_weight = 0;

	min_weight = (*matrix).front();
	max_weight = (*matrix).front();

	// Find min and max weights

	for(int s=0; s<stacks; s++)
	{
		// Horizontal edges
		for(int i=(s*numel_stack); i<= (s*numel_stack)+numel_stack-rows-1; i++)
		{
			// Weighted by matrix
			double weight = abs((*matrix)[i] - (*matrix)[i+rows]);

			if(weight < min_weight )
				min_weight = weight;

			if(weight > max_weight)
				max_weight = weight;
		}

		// Vertical edges
		for(int i=(s*numel_stack); i<= (s*numel_stack)+numel_stack-rows; i+=rows)
		{
			for(int j=0; j<=rows-2; j++)
			{
				double weight = abs((*matrix)[i+j] - (*matrix)[i+j+1]) ;

				if(weight < min_weight )
					min_weight = weight;

				if(weight > max_weight)
					max_weight = weight;
			}
		}
	}

	// Stack edges

	for(int i=0; i<numel-numel_stack; i++)
	{
		// Weighted by matrix
		double weight = abs((*matrix)[i] - (*matrix)[i+numel_stack]);

		if(weight < min_weight )
			min_weight = weight;

		if(weight > max_weight)
			max_weight = weight;
	}

	// Normalization values
	double diff = (max_weight - min_weight);
	// Prevent division by zero
	const double epsilon_diff = 1.0e-17;
	if(diff<epsilon_diff)
		diff = epsilon_diff;

	// Reserve 7 non-zero entries for each column.
	// This is the case for the 6-connected lattice.
	// Always reserve more than less entries (this is very crucial for performance)

	L.reserve(VectorXi::Constant(numel,7));

	//--- Find min and max weights

	for(int s=0; s<stacks; s++)
	{
		// Horizontal edges
		for(int i=(s*numel_stack); i<= (s*numel_stack)+numel_stack-rows-1; i++)
		{
			// .coeffRef save
			// .insert assumes entry is zero, faster but if entry non-zero crash
			// Alternative but a lot slower L.coeffRef(i,i+rows) = 1.0; 	L.coeffRef(i+rows,i) = 1.0;

			// Weighted by matrix
			double weight = abs((*matrix)[i] - (*matrix)[i+rows]);

			weight = (weight - min_weight) / diff;

			weight = exp(-beta*weight) + epsilon; // Compute Gaussian weighting

			L.insert(i,i+rows) = -weight;
			L.insert(i+rows,i) = -weight;
		}

		// Vertical edges
		for(int i=(s*numel_stack); i<= (s*numel_stack)+numel_stack-rows; i+=rows)
		{
			for(int j=0; j<=rows-2; j++)
			{
				double weight = abs((*matrix)[i+j] - (*matrix)[i+j+1]) ;

				weight = (weight - min_weight) / diff;

				weight = exp(-beta*weight) + epsilon; // Compute Gaussian weighting

				L.insert(i+j,i+j+1) = -weight;
				L.insert(i+j+1,i+j) = -weight;
			}
		}
	}

	// Stack edges

	for(int i=0; i<numel-numel_stack; i++)
	{
		// Weighted by matrix
		double weight = abs((*matrix)[i] - (*matrix)[i+numel_stack]);

		weight = (weight - min_weight) / diff;

		weight = exp(-beta*weight) + epsilon; // Compute Gaussian weighting

		L.insert(i,i+numel_stack) = -weight;
		L.insert(i+numel_stack,i) = -weight;
	}
}