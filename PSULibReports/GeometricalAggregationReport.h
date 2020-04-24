#pragma once
#include <string>
#include <vector>

using namespace std;

class GeometricalAggregationReport
{
public:
	void printReport(string datadir);
	void printReport(vector<string> datafiles, string outdir);
private:
	//vector<string> getFilesWithinFolder(string folder);
};

