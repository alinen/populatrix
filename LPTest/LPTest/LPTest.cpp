// LPTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "Eigen/Dense"
#include "lp_lib.h"

void SolveForK(const Eigen::MatrixXd& E, const Eigen::MatrixXd& xBar, Eigen::MatrixXd& k)
{
	lprec *lp;

	int numStates = (int) xBar.rows();
	int numVars = numStates * numStates + numStates; 
	// numVars is a row major matrix

	lp = make_lp(0, numVars); // 0 constraint, 2 variables

	// set objective function
	// create temp variable xic for each state
	double *obj = new double[numVars+1];
	memset(obj, 0, sizeof(double)*(numVars + 1));
	for (int i = 0; i < numStates; i++)
	{
		obj[numStates*numStates + 1 + i] = 1;
	}

	if (!set_obj_fn(lp, obj))
		std::cerr << "Cannot add objective "<< std::endl;

	// add constraints - transitions should be close to 0.5
	for (int i = 0; i < numStates; i++)
	{
		double* row = new double[numVars + 1];
		memset(row, 0, sizeof(double)*(numVars + 1));
		int idx = i * numStates + i;
		int idxc = numStates * numStates + i;

		row[idx + 1] = 1;
		row[idxc + 1] = -1;
		if (!add_constraint(lp, row, LE, 0.5))
			std::cerr << "Cannot add constraint (sum to one) "<< i << std::endl;

		row[idx + 1] = 1;
		row[idxc + 1] = 1;
		if (!add_constraint(lp, row, GE, 0.5))
			std::cerr << "Cannot add constraint (sum to one) "<< i << std::endl;

		delete[] row;
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

	// add constraints - rates are zero when there is no edge

	//double row2[] = { 0, -1, 1 };
	//if (!add_constraint(lp, row2, EQ, 2))
	//	std::cerr << "Cannot add constraint "<< std::endl;

	print_lp(lp);
	printf("%d", solve(lp));
	print_objective(lp);
	print_solution(lp,1);
	print_constraints(lp,1);

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
		
	delete[] obj;
	delete_lp(lp);
}

int main()
{
	Eigen::MatrixXd E(2, 2);
	Eigen::MatrixXd x(2, 1);
	Eigen::MatrixXd k(2, 2);

	E(0, 0) = E(1, 0) = E(0, 1) = E(1, 1) = 1;
	x(0, 0) = 0.9;
	x(1, 0) = 0.1;

	SolveForK(E, x, k);
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
