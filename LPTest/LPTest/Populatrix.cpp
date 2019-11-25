#include "Populatrix.h"
#include <iostream>
#include <fstream>
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

const Eigen::MatrixXd& Populatrix::getDesiredDistribution() const
{
	return _xd;
}

const Eigen::MatrixXi& Populatrix::getEdgeMatrix() const
{
	return _E;
}

const Eigen::MatrixXi& Populatrix::getDurations() const
{
	return _durations;
}

const Eigen::MatrixXd& Populatrix::getRates() const
{
	return _k;
}


void Populatrix::sanityCheck()
{
	assert(_E.rows() == _E.cols()); 
	assert(_xd.size() == _E.rows());
	assert(_durations.size() == 0 || (_durations.size() == _E.rows()));
}

const Eigen::MatrixXd& Populatrix::computeRates()
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
	_k.resize(_E.rows(), _E.cols());
	for (int nodeid = 0; nodeid < n; nodeid++)
	{
		int idx = nodeid;
		if (_durations(nodeid, 0) > 1)
		{
			idx = n + sumDurations + _durations(nodeid, 0) - 1 - 1;
		}
		for (int j = 0; j < n; j++)
		{
			_k(j, nodeid) = ek(j, idx);
		}
		sumDurations += (_durations(nodeid, 0) - 1);
	}
	std::cout << " results\n";
	std::cout << _k << std::endl;
	return _k;
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

void Populatrix::loadRates(const std::string& openName)
{
	int n = 0;

	std::ifstream dFile(openName + "-durations.csv");
	if (!dFile.is_open())
	{
		std::cerr << "Error: cannot find durations for " + openName << std::endl;
		return;
	}
	dFile >> n;
	_durations.conservativeResize(n, 1);
	for (int i = 0; i < n; i++)
	{
		dFile >> _durations(i, 0);
	}
	dFile.close();

	std::ifstream EFile(openName + "-edges.csv");
	if (!EFile.is_open())
	{
		std::cerr << "Error: cannot find edges for " + openName << std::endl;
		return;
	}
	EFile >> n;
	assert(n == _durations.rows());
	_E.conservativeResize(n, n);
	for (int i = 0; i < _E.rows(); i++)
	{
		for (int j = 0; j < _E.cols(); j++)
		{
			EFile >> _E(i, j); 
		}
	}
	EFile.close();

	std::ifstream kFile(openName + "-rates.csv");
	if (!kFile.is_open())
	{
		std::cout << "Cannot find rates for " + openName << std::endl;
		return;
	}
	kFile >> n;
	assert(n == _durations.rows());
	_k.conservativeResize(n, n);
	_xd.conservativeResize(n, 1);
	for (int i = 0; i < _k.rows(); i++)
	{
		for (int j = 0; j < _k.cols(); j++)
		{
			kFile >> _k(i, j);
		}
	}
	for (int i = 0; i < _xd.rows(); i++)
	{
		kFile >> _xd(i, 0);
	}
	kFile.close();
}

void Populatrix::saveRates(const std::string& saveName)
{
	std::ofstream dFile(saveName + "-durations.csv");
	dFile << _durations.size() << std::endl;
	for (int i = 0; i < _durations.size(); i++)
	{
		dFile << _durations(i,0) << std::endl;
	}
	dFile.close();

	std::ofstream EFile(saveName + "-edges.csv");
	EFile << _E.rows() << std::endl;
	for (int i = 0; i < _E.rows(); i++)
	{
		for (int j = 0; j < _E.cols(); j++)
		{
			EFile << _E(i, j) << " "; 
		}
		EFile << std::endl;
	}
	EFile.close();

	std::ofstream kFile(saveName + "-rates.csv");
	kFile << _k.rows() << std::endl;
	for (int i = 0; i < _k.rows(); i++)
	{
		for (int j = 0; j < _k.cols(); j++)
		{
			kFile << _k(i, j) << " "; 
		}
		kFile << std::endl;
	}
	for (int i = 0; i < _xd.size(); i++)
	{
		kFile << _xd(i, 0) << std::endl;
	}
	kFile.close();
}

void Populatrix::loadModel(const std::string& openName)
{
	std::ifstream dFile(openName);
	if (!dFile.is_open())
	{
		std::cerr << "Error: cannot open model: " + openName << std::endl;
		return;
	}
	_keys.clear();
	_areas.clear();
	_activities.clear();
	_sites.clear();

	std::string line;
	while (std::getline(dFile, line))
	{
		std::istringstream linestream(line);
		std::string type;
		linestream >> type;
		if (type == "key")
		{
			std::string name;
			std::string time;
			std::string distStr;
			linestream >> name;
			linestream >> time;
			linestream >> distStr;

			// note: good tokenize appraoch
			std::vector<std::pair<std::string, double>> distribution;
			std::string token;
			std::istringstream distStream(distStr);
			while (std::getline(distStream, token, ','))
			{
				size_t splitPos = token.find(':', 0);
				if (splitPos == std::string::npos)
				{
					std::cerr << "Error parsing key distribution\n";
				}
			    std::string aname = token.substr(0, splitPos);
			    double afraction = stod((token.substr(splitPos+1, token.length() - splitPos)));
				distribution.push_back(std::pair<std::string, double>(aname, afraction));
			}
			unsigned int id = (unsigned int) _keys.size();
			Key key{id, name, time, distribution};
			_keys[id] = key;
		}
		else if (type == "area")
		{
			std::string name;
			std::string locationType;
			linestream >> name;
			linestream >> locationType;

			unsigned int id = (unsigned int) _areas.size();
			Area area{id, name, locationType};
			_areas[id] = area;
		}
		else if (type == "activity")
		{
			std::string name;
			std::string activityType;
			std::string locationType;
			int duration;
			linestream >> name;
			linestream >> activityType;
			linestream >> duration;
			linestream >> locationType;

			unsigned int id = (unsigned int) _activities.size();
			Activity activity{id, name, activityType, locationType, duration};
			_activities[id] = activity;
		}
	}
	initRateModel();
}

void Populatrix::initRateModel() // create rate model from logical model
{
	for (auto item : _areas)
	{
		std::cout << item.second.name << std::endl;
	}

	for (auto item : _activities)
	{
		std::cout << item.second.name << std::endl;
	}

	for (auto item : _keys)
	{
		std::cout << item.second.name << std::endl;
	}
}
