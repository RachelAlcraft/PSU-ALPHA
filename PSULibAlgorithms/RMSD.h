#pragma once

#include<PDBFile.h>
#include<FASTAFile.h>
#include<string>

using namespace std;

class CAlphaPair
{
public:
	Atom* a1;
	Atom* a2;
	CAlphaPair(Atom* aa, Atom* ab)
	{
		a1 = aa;
		a2 = ab;
	}
};


class RMSD
{
public:
	PDBFile* PDB1;
	PDBFile* PDB2;
	FastaFile* Fasta;
	bool Optimise;
	bool Alignment;

private:
	vector<CAlphaPair> _calphaPairs;
public:
	RMSD();
	RMSD(PDBFile* pdb1, PDBFile* pdb2, FastaFile* fasta, bool alignment, bool optimise);
	void SetupCAlphaPairs();
	double calculateRMSD();
	string getAtomMatches();
};


