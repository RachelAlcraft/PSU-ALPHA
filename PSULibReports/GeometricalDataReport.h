#pragma once
#include <PDBFile.h>
class GeometricalDataReport
{
public:
	void printReport(PDBFile* pdb, string filename1, string filename2, string directory, string geodir);
private:
	void printOneReport(PDBFile* pdb, string filename1, string filename2);
};

