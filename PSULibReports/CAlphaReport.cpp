#include "CAlphaReport.h"
#include <LogFile.h>
#include <ProteinManager.h>

void CAlphaReport::printReport(PDBFile* pdb, string fileName)
{//produce data frame report for R reporting
	
	LogFile::getInstance()->writeMessage("Starting C Aplpha distance map");
	
	stringstream report;
	report << "Amino1,Id1,Hydro1,Donicity1,Chemical1,Polar1,Amino2,Id2,Hydro2,Donicity2,Chemical2,Polar2,Distance\n";
	map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdb->pdbCode);
	for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
	{
		Chain* chain = iter->second;
		string chainid = iter->first;
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, chainid);
		vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb->pdbCode, chainid);
		
		for (unsigned int i = 0; i < calphas.size(); ++i)
		{
			Atom* a = calphas[i];
			AminoAcid* aa = aminos[a->aminoId];
			stringstream ss;
			ss << "Distances for: chain " << chain->chainId << ":" << a->atomId << ":" << a->aminoId << ":" << aa->aminoCode;
			LogFile::getInstance()->writeMessage(ss.str());
			
			for (map<string, Chain*>::iterator iterb = iter/* don't double up, start from where we are*/; iterb != chains.end(); ++iterb)
			{
				Chain* chainb = iterb->second;
				vector<Atom*> atomsb = chainb->getCAlphas();
				map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, chainb->chainId);
				unsigned int start = 0;
				bool chainself = false;
				if (iterb == iter)
				{
					start = i;
					chainself = true;
				}
				for (unsigned int j = start; j < atomsb.size(); ++j)
				{
					Atom* b = atomsb[j];
					AminoAcid* ab = aminosb[b->aminoId];
					double distance = a->atomicDistance(b);
					if (distance < 25) //TODO this should be configurable
					{
						report << a->aminoCode << "," << aa->aminoId << ",";
						report << aa->Hydro << ",";
						report << aa->Donicity << ",";
						report << aa->Chemical << ",";
						report << aa->Polar << ",";
						report << b->aminoCode << "," << ab->aminoId << ",";
						report << ab->Hydro << ",";
						report << ab->Donicity << ",";
						report << ab->Chemical << ",";
						report << ab->Polar << ",";
						report << distance << "\n";
						if (!(chainself & (i == j)))
						{
							//and reverse to halve the time
							report << b->aminoCode << "," << ab->aminoId << ",";
							report << ab->Hydro << ",";
							report << ab->Donicity << ",";
							report << ab->Chemical << ",";
							report << ab->Polar << ",";
							report << a->aminoCode << "," << aa->aminoId << ",";
							report << aa->Hydro << ",";
							report << aa->Donicity << ",";
							report << aa->Chemical << ",";
							report << aa->Polar << ",";
							report << distance << "\n";
						}
					}
				}
			}
		}
	}
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}
}

