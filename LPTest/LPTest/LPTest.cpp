// LPTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "time.h"
#include <iostream>
#include <functional>
#include <random>
#include "Eigen/Dense"
#include <unsupported/Eigen/MatrixFunctions>
#include "lp_lib.h"
#include "Populatrix.h"

void SolveForK(const Eigen::MatrixXd& E, const Eigen::MatrixXd& xBar, Eigen::MatrixXd& k)
{
	lprec *lp;

	int numStates = (int) xBar.rows();
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
	double *obj = new double[numVars+1];
	memset(obj, 0, sizeof(double)*(numVars + 1));
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			if (nonZeros(j, 0) > 1)
			{
				obj[numStates*numStates + i*numStates + j + 1] = 1;
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
			int idxc = numStates * numStates + i*numStates + j;

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
			std::cerr << "Cannot add constraint (sum to one) "<< i << std::endl;
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
			row[idx] = xBar(j,0);
		}

		if (!add_constraint(lp, row, EQ, xBar(i,0)))
			std::cerr << "Cannot add constraint (sum to one) "<< i << std::endl;
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
				row[idx] = -xBar(i,0);
			}
		}

		for (int j = 0; j < numStates; j++)
		{
			if (i != j)
			{
				int idx = i * numStates + j + 1;
				row[idx] = xBar(j,0);
			}
		}

		if (!add_constraint(lp, row, EQ, 0))
			std::cerr << "Cannot add constraint (sum to one) "<< i << std::endl;
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
			int idx = i*numStates + j;
			k(i, j) = values[idx];
		}
	}
	delete[] values;
	std::cout << k << std::endl;
	char name[32] = "model.lp";
	write_lp(lp, name);
		
	delete_lp(lp);
}

bool TestK(const Eigen::MatrixXd& k, 
	const Eigen::MatrixXd& xBar, bool verbose = false)
{
	Eigen::MatrixXd rates = k.pow(100);
	if (verbose)
	{
		std::cout << rates << std::endl;
	}

	for (int i = 0; i < rates.rows(); i++)
	{
		for (int j = 0; j < rates.cols(); j++)
		{
			if (std::abs(rates(i, j) - xBar(j, 0)) > 0.000000001)
			{
				return false;
			}
		}
	}

	return true;
}

std::default_random_engine generator((unsigned int) time(0)); 
std::uniform_real_distribution<double> distribution(0, 1); 
auto dice = std::bind(distribution, generator);

double SampleExp(double lambda)
{
	return -log(dice()) / lambda;
}

int ChooseActivity(const Eigen::MatrixXd& k, int current)
{
	auto col = k.col(current);
	double rval = dice();
	double sum = 0;
	for (int i = 0; i < col.size() - 1; i++)
	{
		double next = sum + col(i);
		if (sum <= rval && rval < next)
		{
			return i;
		}
		sum += col(i);
	}
	return (int) col.size() - 1;
}


void TestSimulation(const Eigen::MatrixXd& k,
	const Eigen::MatrixXd& xBar, 
	const Eigen::MatrixXi& durations, 
	int numSteps = 100, 
	int numAgents = 100, 
	bool verbose = false)
{
	// TODO
	struct Agent 
	{
		double time = -1; 
		int activityid = 0;
	};
	Agent* agents = new Agent[numAgents];
	int* counts = new int[xBar.size()];

	double sum = 0;
	Eigen::MatrixXi numStart(xBar.rows(), xBar.cols());
	for (int i = 0; i < numStart.size(); i++)
	{
		sum += xBar(i);
		numStart(i) = (int) (sum * numAgents);
	}
	int activityid = 0;
	for (int i = 0; i < numAgents; i++)
	{
		if (i > numStart(activityid))
		{
			activityid++;
		}
		int time = durations(activityid);
		agents[i].activityid = activityid;
		agents[i].time = time;
	}

	for (int step = 0; step < numSteps; step++)
	{
		for (int i = 0; i < xBar.size(); i++)
		{
			counts[i] = 0;
		}

		for (int i = 0; i < numAgents; i++)
		{
			if (agents[i].time < 0)
			{
				int currentA = agents[i].activityid;
				int a = (int) ChooseActivity(k, currentA);
				int time = durations(a);
				agents[i].activityid = a;
				agents[i].time = time;
			}
			else
			{
				agents[i].time -= 1;
			}
			counts[agents[i].activityid]++;
		}

		std::cout << step << " " ;
		for (int i = 0; i < xBar.size(); i++)
		{
			double percent = counts[i] / ((float)numAgents);
			std::cout << i << ": " << percent << " ";
		}
		std::cout << std::endl;
	}

	delete[] counts;
	delete[] agents;
}

bool TestConvergence(const Eigen::MatrixXd& k,
	const Eigen::MatrixXd& xBar, bool verbose = false)
{
	int n = (int) xBar.rows();
	Eigen::MatrixXd x = Eigen::MatrixXd::Ones(n, 1);

	for (int i = 0; i < n; i++)
	{
		x(i, 0) = 1. / n;
	}

	int steps = 0;
	int maxSteps = 100;
	for (int i = 0; i < maxSteps; i++)
	{
		x = k * x;
		std::cout << x.transpose() << std::endl;
		if ((x - xBar).norm() < 0.01)
		{
			std::cout << "Converged in " << steps << " steps\n";
			break;
		}
		steps++;
	}

	return (steps < maxSteps);
}

// Simple 2 state graph: 0.9 in n1 and 0.1 in n2
void Test1()
{
	Eigen::MatrixXd E(2, 2);
	Eigen::MatrixXd x(2, 1);
	Eigen::MatrixXd k(2, 2);

	E(0, 0) = E(1, 0) = E(0, 1) = E(1, 1) = 1;
	x(0, 0) = 0.9;
	x(1, 0) = 0.1;

	SolveForK(E, x, k);
	TestK(k, x, true);
	TestConvergence(k, x);
}

void Test2()
{
	Eigen::MatrixXd E(3, 3);
	Eigen::MatrixXd x(3, 1);
	Eigen::MatrixXd k(3, 3);

	E(0, 0) = 1; E(0, 1) = 1; E(0, 2) = 1; 
	E(1, 0) = 0; E(1, 1) = 0; E(1, 2) = 1; 
	E(2, 0) = 1; E(2, 1) = 0; E(2, 2) = 1; 

	x(0, 0) = 0.3;
	x(1, 0) = 0.3;
	x(2, 0) = 0.4;

	SolveForK(E, x, k);
	TestK(k, x, true);
	TestConvergence(k, x, true);
}

void Test3()
{
	Eigen::MatrixXd E(4, 4);
	Eigen::MatrixXd x(4, 1);
	Eigen::MatrixXd k(4, 4);

	E(0, 0) = 1; E(0, 1) = 1; E(0, 2) = 0; E(0, 3) = 1;
	E(1, 0) = 1; E(1, 1) = 1; E(1, 2) = 1; E(1, 3) = 0;
	E(2, 0) = 0; E(2, 1) = 1; E(2, 2) = 1; E(2, 3) = 1;
	E(3, 0) = 1; E(3, 1) = 0; E(3, 2) = 1; E(3, 3) = 1;

	x(0, 0) = 0.3;
	x(1, 0) = 0.2;
	x(2, 0) = 0.4;
	x(3, 0) = 0.1;

	SolveForK(E, x, k);
	TestK(k, x, true);
	TestConvergence(k, x, true);

}

void Test4()
{
	Eigen::MatrixXd E(10, 10);
	Eigen::MatrixXd x(10, 1);
	Eigen::MatrixXd k(10, 10);

	// columns are out-going edges
	E << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0;

	x << 0.85,
		0.01,
		0.01,
		0.01,
		0.07,
		0.01,
		0.01,
		0.01,
		0.01,
		0.01;

	SolveForK(E, x, k);
	TestK(k, x, true);
	TestConvergence(k, x, true);
}

void Test5()
{
	Eigen::MatrixXd E(3, 3);
	Eigen::MatrixXd x(3, 1);
	Eigen::MatrixXd k(3, 3);

	E(0, 0) = 1; E(0, 1) = 1; E(0, 2) = 0; 
	E(1, 0) = 1; E(1, 1) = 0; E(1, 2) = 1; 
	E(2, 0) = 1; E(2, 1) = 0; E(2, 2) = 0; 

	x(0, 0) = 0.3;
	x(1, 0) = 0.5;
	x(2, 0) = 0.2;

	SolveForK(E, x, k);
	TestK(k, x, true);
	TestConvergence(k, x, true);
}

void Test6()
{
	Eigen::MatrixXi E(4, 4);
	Eigen::MatrixXd x(4, 1);
	Eigen::MatrixXd k(4, 4);

	E << 1, 1, 0, 1,
	     1, 1, 1, 0,
	     0, 1, 1, 1,
	     1, 0, 1, 1;

	x << 0.3,
	     0.2,
	     0.4,
	     0.1;
	
	Populatrix calculator;
	calculator.setDesiredDistribution(x);
	calculator.setEdgeMatrix(E);
	k = calculator.computeRates();

	TestK(k, x, true);
	TestConvergence(k, x, true);
}

void Test7()
{
	Eigen::MatrixXi E(4, 4);
	Eigen::MatrixXi d(4, 1);
	Eigen::MatrixXd x(4, 1);
	Eigen::MatrixXd k(4, 4);

	E << 1, 0, 0, 1, 
		 1, 0, 0, 0, 
		 0, 1, 1, 0, 
		 0, 0, 1, 0; 

	x << 0.85,
		 0.03,
		 0.07,
		 0.05;

	d << 1,
		3,
		1,
		5;

	Populatrix calculator;
	calculator.setDesiredDistribution(x);
	calculator.setDurations(d);
	calculator.setEdgeMatrix(E);
	k = calculator.computeRates();
	calculator.saveRates("test");

	TestK(k, x, true);
	TestConvergence(k, x, true);
	TestSimulation(k, x, d, 1000, 100, true);
}

void Test8()
{
	Populatrix calculator;
	calculator.loadRates("test");

	Eigen::MatrixXi E = calculator.getEdgeMatrix();
	Eigen::MatrixXi d = calculator.getDurations();
	Eigen::MatrixXd x = calculator.getDesiredDistribution();
	Eigen::MatrixXd k = calculator.getRates();

	std::cout << E << std::endl;
	std::cout << d << std::endl;
	std::cout << x << std::endl;
	std::cout << k << std::endl;

	TestK(k, x, true);
	TestConvergence(k, x, true);
	TestSimulation(k, x, d, 1000, 100, true);
}

int main()
{
	Populatrix p;
	p.loadModel("model.csv");
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
