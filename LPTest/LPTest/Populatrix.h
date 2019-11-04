#pragma once
#include <Eigen/Dense>
#include <vector>

class Populatrix
{
public:
	Populatrix();
	virtual ~Populatrix();
	
	void setDesiredDistribution(const Eigen::MatrixXd& xd);
	void setEdgeMatrix(const Eigen::MatrixXi& E);
	void setDurations(const Eigen::MatrixXi& durations);
	void setTravelTimes(const Eigen::MatrixXi& D);

	void computeRates(Eigen::MatrixXd& k);

	/*
	int getNumActivities() const;
	std::string getName(int id) const;
	float getDuration(int id) const;
	int getId(const std::string& name) const;
	std::vector<int> getNextActivities() const;
	std::vector<int> getPrevActivities() const;
	float getTransitionRate(int id) const;
	*/

protected:

	void expandGraph(Eigen::MatrixXi& E, Eigen::MatrixXd& xd);
	void computeExpandedRates(const Eigen::MatrixXi& E, 
		const Eigen::MatrixXd& xd, Eigen::MatrixXd& k);
	void sanityCheck();

	struct Activity
	{
		int id;
		std::string name;
		int duration;
	};
	Eigen::MatrixXi _durations;
	Eigen::MatrixXi _E; // edge matrix
	Eigen::MatrixXd _xd; // desired activity distribution
	Eigen::MatrixXd _k; // rates matrix
	Eigen::MatrixXi _D; // travel time matrix
	// todo: need areas where activities can be performed
};
