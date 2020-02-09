#pragma once

#include <string>
#include <RMSD.h>


using namespace std;


class RMSDReport
{
public:
	void printReport(RMSD* rmsd, string fileName, bool optimised, string fileroot);
};

