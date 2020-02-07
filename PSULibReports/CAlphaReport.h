#pragma once
#include <PDBFile.h>
#include <string>

using namespace std;

class CAlphaReport
{
public:
	void printReport(PDBFile* pdb, string fileName);
};

