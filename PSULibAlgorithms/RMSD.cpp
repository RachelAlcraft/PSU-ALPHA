#include "RMSD.h"
#include <LogFile.h>
#include <ProteinManager.h>
#include <algorithm>
#include <GeoShape.h>

RMSD::RMSD()
{
	PDB1 = nullptr;
	PDB2 = nullptr;
	Fasta = nullptr;
	Optimise = false;
	Alignment = false;
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
			_geo1.addCoords(calphas1[i]->shifted_coords);
			_geo1.addCoords(calphas2[i]->shifted_coords);
			LogFile::getInstance()->writeMessage("RMSD Match: " + calphas1[i]->getDescription() + " " + calphas2[i]->getDescription());
		}
	}
	else
	{//use alignment file to match off

	}
}

double RMSD::calculateOptimalRMSD()
{
	//First calculation is after the most obvious spatial alignment		
	GeoTripod tri1 = _geo1.getTripod(1,1); // 1 is the best solution
	GeoTripod tri2 = _geo2.getTripod(1,1);// 1 is the best solution
	GeoTransformation gt = tri1.getTransformation(tri2);
	PDB2->applyTransformation(gt);	
	
	return calculateRMSD();
	
	//TODO later I can apply more shifts and see if this is not optimal
	// Do this by iterating through the tripod solutions
}

double RMSD::calculateRMSD() // this may be iteratively called from an optimise function and needs to be fast
{			
	double rmsd = 0.0;
	for (unsigned int i = 0; i < _calphaPairs.size(); ++i)
	{//We are always using the shifted coordinates for this.
		GeoCoords aCoords = _calphaPairs[i].a1->shifted_coords;
		GeoCoords bCoords = _calphaPairs[i].a2->shifted_coords;
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
	for (unsigned int i = 0; i < _calphaPairs.size(); ++i)	
		ss << _calphaPairs[i].a1->getDescription() << ":" << _calphaPairs[i].a2->getDescription() << "\n";;
	return ss.str();
}

