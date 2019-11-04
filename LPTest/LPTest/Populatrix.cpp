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

void Populatrix::setDurations(const Eigen::MatrixXi& durations)
{
	_durations = durations;
}

void Populatrix::sanityCheck()
{
	assert(_E.rows() == _E.cols()); 
	assert(_xd.size() == _E.rows());
	assert(_durations.size() == 0 || (_durations.size() == _E.rows()));
}

void Populatrix::computeRates(Eigen::MatrixXd& k)
{
	sanityCheck();

	Eigen::MatrixXi eE;  // expanded E
	Eigen::MatrixXd eXd; // expanded xd
	expandGraph(eE, eXd);

	Eigen::MatrixXd ek; // expanded k
	computeExpandedRates(eE, eXd, ek);

	// copy over rates to simple k
	// the rates will be contained in the last node of each chain
	int n = (int) _E.rows();
	int sumDurations = 0;
	k.resize(_E.rows(), _E.cols());
	for (int nodeid = 0; nodeid < n; nodeid++)
	{
		int idx = nodeid;
		if (_durations(nodeid, 0) > 1)
		{
			idx = n + sumDurations + _durations(nodeid, 0) - 1 - 1;
		}
		for (int j = 0; j < n; j++)
		{
			k(j, nodeid) = ek(j, idx);
		}
		sumDurations += (_durations(nodeid, 0) - 1);
	}
	std::cout << " results\n";
	std::cout << k << std::endl;
}

static void RemoveEdge(Eigen::MatrixXi& M, int i, int j)
{
	M(i, j) = 0;
}

static void AddEdge(Eigen::MatrixXi& M, int i, int j)
{
	M(i, j) = 1;
}

static int AddNode(Eigen::MatrixXi& M)
{
	int n = (int) M.rows();
	M.conservativeResize(n + 1, n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		M(n, i) = 0;
		M(i, n) = 0;
	}
	return n;
}

void Populatrix::expandGraph(
	Eigen::MatrixXi& E, Eigen::MatrixXd& xd)
{
	E = _E;
	xd = _xd;
	std::cout << "TEST 1 \n" << E << std::endl;
	std::cout << "TEST 1 \n" << xd << std::endl;

	if (_durations.size() > 0) 
	{
		int n = (int) _durations.rows();
		for (int nodeid = 0; nodeid < n; nodeid++) 
		{
			int d = _durations(nodeid, 0);
			if (d > 1)
			{
				// save existing edges
				auto col = _E.col(nodeid);

				// set these edges to zero
				for (int row = 0; row < n; row++)
				{
					RemoveEdge(E, row, nodeid);
				}

				// create chain
				// extend xd too
				xd.conservativeResize(xd.rows()+d-1,1);
				xd(nodeid) = _xd(nodeid) / d;
				int previd = nodeid;
				int id = 0;
				for (int row = 0; row < d-1; row++)
				{
					id = AddNode(E);
					AddEdge(E, id, previd);

					xd(id) = _xd(nodeid) / d;
					previd = id;
					std::cout << "LOOP " << row << " \n" << E << std::endl;
					std::cout << "LOOP " << row << " \n" << xd << std::endl;
				}

				// last node has same edges as original node
				for (int row = 0; row < col.size(); row++)
				{
					if (col(row) == 1)
					{
						AddEdge(E, row, id);
					}
				}
				std::cout << "TEST 2 \n" << E << std::endl;
				std::cout << "TEST 2 \n" << xd << std::endl;

			}
		}
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

