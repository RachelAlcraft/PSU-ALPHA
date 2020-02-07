#include "RMSD.h"
#include <LogFile.h>
#include <ProteinManager.h>
#include <algorithm>

RMSD::RMSD()
{
	PDB1 = nullptr;
	PDB2 = nullptr;
	Fasta = nullptr;
	Optimise = false;
}
RMSD::RMSD(PDBFile* pdb1, PDBFile* pdb2, FastaFile* fasta, bool alignment, bool optimise)
{
	PDB1 = pdb1;
	PDB2 = pdb2;
	Fasta = fasta;
	Alignment = alignment;
	Optimise = optimise;
	SetupCAlphaPairs();
}

void RMSD::SetupCAlphaPairs()
{//For a very first implementation this will just be the best matching cAlphas
	vector<Atom*> calphas1 = ProteinManager::getInstance()->getCAlphas(PDB1->pdbCode);
	vector<Atom*> calphas2 = ProteinManager::getInstance()->getCAlphas(PDB2->pdbCode);
	
	int maxposs_calphas = min(calphas1.size(), calphas2.size());

	if (!Alignment)
	{
		for (int i = 0; i < maxposs_calphas; ++i)
		{			
			_calphaPairs.push_back(CAlphaPair(calphas1[i],calphas2[i]));
			LogFile::getInstance()->writeMessage("RMSD Match: " + calphas1[i]->getDescription() + " " + calphas2[i]->getDescription());
		}
	}
	else
	{//use alignment file to match off

	}
}

double RMSD::calculateRMSD() // this may be iteratively called from an optimise function and needs to be fast
{			
	double rmsd = 0.0;
	for (unsigned int i = 0; i < _calphaPairs.size(); ++i)
	{//We are always using the shifted coordinates for this.
		Coordinates aCoords = _calphaPairs[i].a1->shifted_coords;
		Coordinates bCoords = _calphaPairs[i].a2->shifted_coords;		
		double distance = 0.0;
		distance += pow((aCoords.x - bCoords.x), 2);
		distance += pow((aCoords.y - bCoords.y), 2);
		distance += pow((aCoords.z - bCoords.z), 2);
		rmsd += distance;
	}	
	return rmsd;
}

string RMSD::getAtomMatches()
{
	stringstream ss;	
	for (int i = 0; i < _calphaPairs.size(); ++i)	
		ss << _calphaPairs[i].a1->getDescription() << ":" << _calphaPairs[i].a2->getDescription() << "\n";;
	return ss.str();
}

