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
			_geo2.addCoords(calphas2[i]->shifted_coords);
			LogFile::getInstance()->writeMessage("RMSD Match: " + calphas1[i]->getDescription() + " " + calphas2[i]->getDescription());
		}
	}
	else
	{//use alignment file to match off

	}
}



string RMSD::calculateRMSD() // this may be iteratively called from an optimise function and needs to be fast
{
	stringstream ss;
	if (Optimise)
	{
		ss << "Optimised report\n";
		ss << "Initial Value=" << calculateOneRMSD() << "\n";
		//Dummy attempt at optimisation TODO ALSO TODO if any ever come back as 0 we can stop!
		double best = 0;
		unsigned int hval = 0;
		unsigned int ival = 0;
		unsigned int jval = 0;
		unsigned int kval = 0;		
		unsigned int MAXITER = 3; //TESTING ONLY DO IT ONCE		
		for (unsigned int h = 1; h < MAXITER; ++h)
		{
			for (unsigned int i = 1; i < MAXITER; ++i)
			{
				for (unsigned int j = 1; j < MAXITER; ++j)
				{
					for (unsigned int k = 1; k < MAXITER; ++k)
					{						
						double rmsd = calculateOptimalRMSD(h, i, j, k);
						ss << "Opt: h=" << h << " i=" << i << " j=" << j << " k=" << k << " rmsd=" << rmsd << "\n";
						if (ival == 0 || rmsd < best)
						{
							hval = h;
							ival = i;
							jval = j;
							kval = k;
							//orientation = orient;
							best = rmsd;
						}						
					}
				}
			}
		}
		
		//So we have done a minimum optimisation and we will take the best, calculating it again TODO because I haven't saved it
		double rmsd = calculateOptimalRMSD(hval,ival, jval,kval);		
		ss << "Optimised report: RMSD Value=" << rmsd;
	}
	else
	{
		ss << "Non optimised report: RMSD Value=" << calculateOneRMSD();
	}
	return ss.str();	
}
double RMSD::calculateOptimalRMSD(int h,int i, int j, int k/*, int orientation*/)
{
	stringstream ss;
	GeoTripod tri1, tri2;
	_geo1.makeTripod(tri1,h, i); // 1 is the best solution
	_geo2.makeTripod(tri2,j, k);// 1 is the best solution
	//For now I am moving both structures onto the orgin to compare them as I have failed to move one on to the other :-( TODO I could just also go backwards but for now this will do
	GeoTransformations* gt1 = tri1.getTransformation(tri1/*,orientation*/);	
	GeoTransformations* gt2 = tri2.getTransformation(tri2/*, orientation*/);
	PDB1->applyTransformation(gt1);
	PDB2->applyTransformation(gt2);
	double val = calculateOneRMSD();		
	return val;
}
double RMSD::calculateOneRMSD() // this may be iteratively called from an optimise function and needs to be fast
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

