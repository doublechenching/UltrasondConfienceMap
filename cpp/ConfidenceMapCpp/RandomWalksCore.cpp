#include "RandomWalksCore.h"
#include "SparseSolverFactory.h"
#include <iostream>

RandomWalksCore::RandomWalksCore():
sparseSolver(0),
isLaplaceAvailable(false)
{
	solverFactory = new SparseSolverFactory();
	sparseSolver = solverFactory->createSolver("Eigen-LLT");
}

RandomWalksCore::~RandomWalksCore()
{
	if(sparseSolver!=0)
		delete sparseSolver;
	sparseSolver = 0;

	if(solverFactory!=0)
		delete solverFactory;
	solverFactory = 0;
}

void RandomWalksCore::setSolver(std::string solver, int iterations, double tolerance)
{
	if(sparseSolver!=0)
		delete sparseSolver;

	sparseSolver = solverFactory->createSolver(solver, iterations, tolerance);
}

std::vector<double> RandomWalksCore::solve()
{
	if(this->matrix->size() <= 0 || this->seeds->size() <= 0 || this->labels->size() <=0)
	{
		std::cout << "ERROR: External inputs (matrix, seeds, labels, etc.) not available\n";
		std::vector<double> nullSolution;
		return nullSolution;
	}

	if(!isLaplaceAvailable)
	{
		//std::cout << "Solving..." << std::endl;
		// Assemble Laplacian
		//std::cout << "Assemble Edges..." << std::endl;
		this->assembleLaplacianEdges(); // First edges, because those are implementation dependent
		//std::cout << "Assemble Degree..." << std::endl;
		this->assembleDegree(); // Then complete with degree
		//std::cout << "Unique indices..." << std::endl;
		this->generateUniqueIndices();
		isLaplaceAvailable = true;
	}
	//std::cout << "Assemble Ax = b..." << std::endl;
	this->assemble_Lu_b();

	//std::cout << "Boundary conditions applied..." << std::endl;
	//std::cout << "Solving sparse Ax=b..." << std::endl;
	return this->sparseSolver->solve_Ax_b(Lu,b,numel,uidx,labels,seeds,active_label);
}

void RandomWalksCore::setLabeling(const std::vector<int>  * seeds, const std::vector<int> * labels, int active_label, int num_labels)
{
	this->seeds = seeds;
	this->labels = labels;
	this->active_label = active_label;
	this->num_labels = num_labels;
	this->isLaplaceAvailable = false;
}

void RandomWalksCore::assemble_Lu_b()
{
	const int n = seeds->size();// # marked nodes
	const int q = numel-seeds->size(); // # unmarked nodes

	// Permutation matrix to get UNMARKED COLUMNS
	SparseMatrix<double> Cu(numel,q);

	for (int i=0; i<q; i++)
	{
		Cu.insert(uidx[i],i) = 1.0; // InsertBack because we know for this column that we place the last entry into the inner vector Eigen storage format. Its more efficient than insert.
	}

	// Cu.transpose() means actually rows unmarked. Once again A^T C A is encountered.
	this->Lu = Cu.transpose() * L * Cu;

	// Permutation matrix to get MARKED COLUMNS
	SparseMatrix<double> Cm(numel,n);

	for (int i=0; i<n; i++)
	{
		Cm.insert((*seeds)[i],i) = 1.0;
	}
	SparseMatrix<double> Bt(q,n);

	Bt = Cu.transpose() * L * Cm;

	SparseMatrix<double> M(n,num_labels);

	for (int i=0; i<num_labels; i++)
	{
		for(int j=0; j<n; j++)
		{
			if((*labels)[j] == i)
				M.insert(j,i) = 1.0;
		}
	}

	SparseVector<double> CL(M.cols());
	CL.insert(active_label) = 1.0;
	this->b = -Bt*M*CL;
}

void RandomWalksCore::generateUniqueIndices()
{
	const int n = seeds->size();// # marked nodes
	const int q = numel-seeds->size(); // # unmarked nodes

	std::vector<int> uidx_temp(numel);

	// Set all node indices
	for(int i=0; i<numel; i++)
	{
		uidx_temp[i] = i;
	}

	// Flag marked nodes with negative sign (last -1 in case of 0th node being a seed
	for(int i=0; i<n; i++)
	{
		uidx_temp[(*seeds)[i]] = -uidx_temp[(*seeds)[i]]-1;
	}

	// Vecotr of unmarked nodes. Marked nodes are in seeds vector
	uidx.reserve(q);

	for(int i=0; i<numel; i++)
	{
		if(uidx_temp[i] >= 0)
			uidx.push_back(uidx_temp[i]);
	}
}

void RandomWalksCore::assembleDegree()
{
	// Filling the degree (diagonal) of Lacency matrix
	// Equals the sum of the row or column. Both are equivalent but looping throught the inner loop
	// is more efficient for Eigen. Therefore, since we declared an column-major matrix we sum up the column entries.

	double degree = 0.0;

	for (int k=0; k<L.outerSize(); ++k)
	{
		degree = 0.0;
		for (SparseMatrix<double>::InnerIterator it(L,k); it; ++it)
		{
			degree += it.value();
		}
		L.insert(k,k) = abs(degree);
	}
}