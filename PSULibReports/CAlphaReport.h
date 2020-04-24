#pragma once
#include <PDBFile.h>
#include <string>

using namespace std;

class CAlphaReport
{
public:
	void printReport(PDBFile* pdb, string occupant, string chain1, string chain2, string fileName);
	void printSingleChainReport(PDBFile* pdb, string occupant, string chain1, string chain2, string fileName);
	void printMultiReport(PDBFile* pdb1, PDBFile* pdb2,string occupant, string fileName, bool shifted);
private:
	string getRow(AminoAcid* aa, Atom* a, AminoAcid* ab, Atom* b, double distance, bool singleChain);
	string getHeader();
};

