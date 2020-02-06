#pragma once

#include <PDBFile.h>
#include <string>

using namespace std;

class RamaReport
{
public:
	void printReport(PDBFile* pdb, string filename);
};

