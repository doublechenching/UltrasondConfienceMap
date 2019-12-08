#include "ConfidenceMaps2D.h"

ConfidenceMaps2D::ConfidenceMaps2D()
{
}

ConfidenceMaps2D::~ConfidenceMaps2D()
{
}

void ConfidenceMaps2D::setMatrix2D(const std::vector<double> * matrix, int rows, int cols, bool normalizeValues)
{
	this->matrix = matrix;
	this->rows = rows;
	this->cols = cols;
	this->numel = rows*cols;

	std::vector<double> * matrix_ptr = const_cast<std::vector<double>* >(matrix);

	if(normalizeValues)
	{
		double min_mat;
		double max_mat;

		min_mat = (*matrix).front();
		max_mat = (*matrix).front();

		//Find min-max
		for(int i=0; i<matrix->size();i++)
		{
			double val = (*matrix)[i];
			if( val< min_mat )
				min_mat = val;

			if(val > max_mat)
				max_mat = val;
		}

		// Normalization values
		double diff = (max_mat - min_mat);
		// Prevent division by zero
		const double epsilon_diff = 1.0e-17;
		if(diff<epsilon_diff)
			diff = epsilon_diff;

		for(int i=0; i<matrix_ptr->size();i++)
			(*matrix_ptr)[i] = ((*matrix_ptr)[i] - min_mat) /diff;
	}

	std::vector<double> dist_penalty;
	dist_penalty.reserve(rows);
	for(int i=0; i<rows; i++)
	{
		double val = 1.0- (exp(- alpha * (double(i)/double(rows))));
		dist_penalty.push_back(val);
	}

	// Apply distance weighting to image
	for(int j=0;j<cols;j++){
		for(int i=1; i<rows; i++){
			 (*matrix_ptr)[i+j*rows] = (*matrix_ptr)[i+j*rows]*dist_penalty[i];
		}
	}
}

void ConfidenceMaps2D::assembleLaplacianEdges()
{
	//const T epsilon = std::numeric_limits<T>::epsilon(); // Numerical limit to avoid division by zero and add to the weights to avoid zero weights in Laplacian

	// Numerical limit to avoid division by zero and add to the weights to avoid zero weights in Laplacian
	const double epsilon = 1.0e-5;

	this->L.resize(numel,numel);

	// Find min and max weights

	double min_weight;
	double max_weight;

	min_weight = (*matrix).front();
	max_weight = (*matrix).front();

	// Horizontal edges
	for(int i=0; i<= numel-rows-1; i++)
	{
		double weight = abs((*matrix)[i] - (*matrix)[i+rows]);
		weight += gamma;

		if(weight < min_weight )
			min_weight = weight;

		if(weight > max_weight)
			max_weight = weight;
	}

	// Vertical edges
	for(int i=0; i<= numel-rows; i+=rows)
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

	// ---Find min and max weights

	// Normalization values
	double diff = (max_weight - min_weight);
	// Prevent division by zero
	const double epsilon_diff = 1.0e-17;
	if(diff<epsilon_diff)
		diff = epsilon_diff;

	// Reserve 5 non-zero entries for each column.
	// This is the case for the 4-connected lattice.
	// Always reserve more than less entries (this is very crucial for performance)
	L.reserve(VectorXi::Constant(numel,5));

	// Horizontal edges
	for(int i=0; i<= numel-rows-1; i++)
	{
		// .coeffRef save
		// .insert assumes entry is zero, faster but if entry non-zero crash
		// Alternative but a lot slower L.coeffRef(i,i+rows) = 1.0; 	L.coeffRef(i+rows,i) = 1.0;

		// Weighted by matrix
		double weight = abs((*matrix)[i] - (*matrix)[i+rows]);
		weight += gamma;
		weight = (weight - min_weight) / diff;
		weight = exp(-beta*weight) + epsilon; // Compute Gaussian weighting

		L.insert(i,i+rows) = -weight;
		L.insert(i+rows,i) = -weight;
	}

	// Vertical edges
	for(int i=0; i<= numel-rows; i+=rows)
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