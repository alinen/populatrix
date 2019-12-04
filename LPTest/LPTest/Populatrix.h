#pragma once
#include <Eigen/Dense>
#include <vector>
#include <map>

class Populatrix
{
public:
	Populatrix();
	virtual ~Populatrix();
	
	void loadModel(const std::string& openName);

	void loadRates(const std::string& openName);
	void saveRates(const std::string& saveName);

	void setDesiredDistribution(const Eigen::MatrixXd& xd);
	void setEdgeMatrix(const Eigen::MatrixXi& E);
	void setDurations(const Eigen::MatrixXi& durations);

	const Eigen::MatrixXd& getDesiredDistribution() const;
	const Eigen::MatrixXi& getEdgeMatrix() const;
	const Eigen::MatrixXi& getDurations() const;
	const Eigen::MatrixXd& getRates() const;

	const Eigen::MatrixXd& computeRates();

protected:

	void expandGraph(Eigen::MatrixXi& E, Eigen::MatrixXd& xd);
	void computeExpandedRates(const Eigen::MatrixXi& E, 
		const Eigen::MatrixXd& xd, Eigen::MatrixXd& k);
	void sanityCheck();
	void initRateModel(); // create rate model from logical model

	// rate model: matrices based on activity,site ids
	Eigen::MatrixXi _durations;
	Eigen::MatrixXi _E; // edge matrix
	Eigen::MatrixXd _xd; // desired activity distribution
	Eigen::MatrixXd _k; // desired activity distribution

	// logical model: human annotation and keys
	struct Site
	{
		unsigned int activityId;
		unsigned int areaId;
		friend bool operator > (const Site& lhs, const Site& rhs)
		{
			return lhs.activityId > rhs.activityId &&
				lhs.areaId >= rhs.areaId;
		}
		friend bool operator < (const Site& lhs, const Site& rhs)
		{
			return !(lhs > rhs); 
		}
	};

	struct Activity
	{
		unsigned int id;      // unique id
		std::string name;     // ex: Dance
		std::vector<std::string> areas;     // area or area cateory, ex. Plaza
		int duration;         // duration in 'ratio' units	
	};

	struct Area
	{
		unsigned int id;
		std::string name;  // e.g. Area0
		std::string type;  // e.g. Plaza
		std::vector<std::string> activities;
	};

	struct Key
	{
		unsigned int id;
		std::string name;
		std::string time;
		std::vector<std::pair<std::string, double>> distribution;
	};

	std::map<int, Key> _keys;
	std::map<int, Activity> _activities;
	std::map<int, Area> _areas;
	std::vector<Site> _sites;
};
