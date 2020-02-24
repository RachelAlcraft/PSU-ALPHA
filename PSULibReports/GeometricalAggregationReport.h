#pragma once
#include <string>
#include <vector>

using namespace std;

class GeometricalAggregationReport
{
public:
	void printReport(string datadir, string filename);
private:
	vector<string> getFilesWithinFolder(string folder);
};

