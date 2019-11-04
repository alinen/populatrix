#include "Populatrix.h"
#include <iostream>
#include "lp_lib.h"

Populatrix::Populatrix()
{
}

Populatrix::~Populatrix()
{
}

void Populatrix::setDesiredDistribution(const Eigen::MatrixXd& xd)
{
	_xd = xd;
}

void Populatrix::setEdgeMatrix(const Eigen::MatrixXi& E)
{
	_E = E;
}

void Populatrix::setTravelTimes(const Eigen::MatrixXi& D)
{
	_D = D;
}

void Populatrix::setDurations(const Eigen::MatrixXi& durations)
{
	_durations = durations;
}

void Populatrix::sanityCheck()
{
	assert(_E.rows() == _E.cols()); 
	assert(_xd.size() == _E.rows());
	assert(_durations.size() == 0 || (_durations.size() == _E.rows()));
	assert(_D.size() == 0 || (_D.rows() == _E.rows()));
}

void Populatrix::computeRates(Eigen::MatrixXd& k)
{
	sanityCheck();

	Eigen::MatrixXi eE;  // expanded E
	Eigen::MatrixXd eXd; // expanded xd
	expandGraph(eE, eXd);

	Eigen::MatrixXd ek; // expanded k
	computeExpandedRates(eE, eXd, ek);

	// todo: copy over rates to simple k
	// this should remove all extra edges used to model timing
	k = ek;
}


void Populatrix::expandGraph(
	Eigen::MatrixXi& E, Eigen::MatrixXd& xd)
{
	if (_D.size() == 0 && _durations.size() == 0)
	{
		E = _E;
		xd = _xd;
		return;
	}

}

void Populatrix::computeExpandedRates(
	const Eigen::MatrixXi& E, const Eigen::MatrixXd& xd, 
	Eigen::MatrixXd& k)
{
	lprec *lp;

	int numStates = (int)xd.rows();
	int numVars = 2 * numStates * numStates;
	// numVars is a row major matrix

	lp = make_lp(0, numVars); // 0 constraint, 2 variables

	// add constraints - rates are zero when there is no edge
	Eigen::MatrixXd nonZeros = Eigen::MatrixXd::Zero(numStates, 1);
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			if (E(i, j) == 0)
			{
				double* row = new double[numVars + 1];
				memset(row, 0, sizeof(double)*(numVars + 1));
				row[i*numStates + j + 1] = 1;
				if (!add_constraint(lp, row, EQ, 0))
					std::cerr << "Cannot add constraint " << std::endl;
				delete[] row;
			}
			else
			{
				nonZeros(j, 0) += 1;
				double* row = new double[numVars + 1];
				memset(row, 0, sizeof(double)*(numVars + 1));
				row[i*numStates + j + 1] = 1;
				if (!add_constraint(lp, row, GE, 0.001))
					std::cerr << "Cannot add constraint " << std::endl;
				delete[] row;
			}
		}
	}

	// set objective function
	// create temp variable xic for each variable
	double *obj = new double[numVars + 1];
	memset(obj, 0, sizeof(double)*(numVars + 1));
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			if (nonZeros(j, 0) > 1)
			{
				obj[numStates*numStates + i * numStates + j + 1] = 1;
			}
		}
	}
	if (!set_obj_fn(lp, obj))
		std::cerr << "Cannot add objective " << std::endl;
	delete[] obj;

	// add constraints - transitions for non-zero rates should be close to 1/n
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			if (E(i, j) == 0) continue;
			if (nonZeros(j, 0) < 2) continue;

			double* row = new double[numVars + 1];
			memset(row, 0, sizeof(double)*(numVars + 1));
			int idx = i * numStates + j;
			int idxc = numStates * numStates + i * numStates + j;

			double preferredRate = 1. / nonZeros(j, 0);

			row[idx + 1] = 1;
			row[idxc + 1] = -1;
			if (!add_constraint(lp, row, LE, preferredRate))
				std::cerr << "Cannot add constraint (sum to one) " << i << std::endl;

			row[idx + 1] = 1;
			row[idxc + 1] = 1;
			if (!add_constraint(lp, row, GE, preferredRate))
				std::cerr << "Cannot add constraint (sum to one) " << i << std::endl;

			delete[] row;
		}
	}

	// add constraints - columns sum to one
	for (int i = 0; i < numStates; i++)
	{
		// column i should sum to one
		double* row = new double[numVars + 1];
		memset(row, 0, sizeof(double)*(numVars + 1));
		for (int j = 0; j < numStates; j++)
		{
			int idx = i + j * numStates + 1;
			row[idx] = 1;
		}

		if (!add_constraint(lp, row, EQ, 1))
			std::cerr << "Cannot add constraint (sum to one) " << i << std::endl;
		delete[] row;
	}

	// add constraints - transitions do not change distribution at steady state
	for (int i = 0; i < numStates; i++)
	{
		double* row = new double[numVars + 1];
		memset(row, 0, sizeof(double)*(numVars + 1));
		for (int j = 0; j < numStates; j++)
		{
			int idx = i * numStates + j + 1;
			row[idx] = xd(j, 0);
		}

		if (!add_constraint(lp, row, EQ, xd(i, 0)))
			std::cerr << "Cannot add constraint (sum to one) " << i << std::endl;
		delete[] row;
	}

	// add constraints - rates at steady state do not change
	for (int i = 0; i < numStates; i++)
	{
		// column i should sum to one
		double* row = new double[numVars + 1];
		memset(row, 0, sizeof(double)*(numVars + 1));
		for (int j = 0; j < numStates; j++)
		{
			if (i != j)
			{
				int idx = i + j * numStates + 1;
				row[idx] = -xd(i, 0);
			}
		}

		for (int j = 0; j < numStates; j++)
		{
			if (i != j)
			{
				int idx = i * numStates + j + 1;
				row[idx] = xd(j, 0);
			}
		}

		if (!add_constraint(lp, row, EQ, 0))
			std::cerr << "Cannot add constraint (sum to one) " << i << std::endl;
		delete[] row;
	}

	//print_lp(lp);
	printf("%d", solve(lp));
	print_objective(lp);
	//print_solution(lp,1);
	//print_constraints(lp,1);

	double* values = new double[numVars];
	get_variables(lp, values);
	k.resize(numStates, numStates);// maybe not necessary
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			int idx = i * numStates + j;
			k(i, j) = values[idx];
		}
	}
	delete[] values;
	std::cout << k << std::endl;
	char name[32] = "model.lp";
	write_lp(lp, name);

	delete_lp(lp);
}

