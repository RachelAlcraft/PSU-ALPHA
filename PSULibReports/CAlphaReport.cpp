#include "CAlphaReport.h"
#include <LogFile.h>
#include <ProteinManager.h>

void CAlphaReport::printReport(PDBFile* pdb, string occupant, string chain1, string chain2, string fileName)
{//produce data frame report for R reporting
	
	LogFile::getInstance()->writeMessage("Starting C Aplpha distance map");
	if (chain1 != "" && chain2 != "")
	{
		printSingleChainReport(pdb, occupant, chain1, chain2, fileName);
	}
	else
	{
		stringstream report;
		//stringstream report_gnuX;
		//stringstream report_gnuY;
		report << getHeader() << "\n";
		map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdb->pdbCode, occupant);
		for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
		{
			Chain* chain = iter->second;
			string chainid = iter->first;
			map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, occupant,chainid);
			vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb->pdbCode, occupant,chainid);

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
					map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, "A",chainb->chainId);
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
						double distance = a->atomicDistance(b,false);
						if (distance < 25) //TODO this should be configurable
						{
							report << getRow(aa, a, ab, b, distance,false) << "\n";							
							if (!(chainself & (i == j)))//and reverse to halve the time														
							{
								report << getRow(ab, b, aa, a, distance, false) << "\n";
								//report_gnuX << ;
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
}

void CAlphaReport::printSingleChainReport(PDBFile* pdb, string occupant,string chain1, string chain2, string fileName)
{//produce data frame report for R reporting

	stringstream report;
	report << getHeader() << "\n";
	Chain* pchain1 = ProteinManager::getInstance()->getOrAddChain(pdb->pdbCode,occupant,chain1);
	Chain* pchain2 = ProteinManager::getInstance()->getOrAddChain(pdb->pdbCode, occupant, chain2);
	map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, occupant, chain1);
	vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb->pdbCode, occupant, chain1);

	for (unsigned int i = 0; i < calphas.size(); ++i)
	{
		Atom* a = calphas[i];
		AminoAcid* aa = aminos[a->aminoId];
		stringstream ss;
		ss << "Distances for: chain " << chain1 << ":" << a->atomId << ":" << a->aminoId << ":" << aa->aminoCode;
		LogFile::getInstance()->writeMessage(ss.str());
				
		vector<Atom*> atomsb = pchain2->getCAlphas();
		map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, occupant, chain2);
		for (unsigned int j = 0; j < atomsb.size(); ++j)
		{
			Atom* b = atomsb[j];
			AminoAcid* ab = aminosb[b->aminoId];
			double distance = a->atomicDistance(b,false);
			if (distance < 25) //TODO this should be configurable
			{				
				report << getRow(aa, a, ab, b, distance,true) << "\n";
			}
		}			
	}
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}
	
}

void CAlphaReport::printMultiReport(PDBFile* pdb1, PDBFile* pdb2, string occupant, string fileName, bool shifted)
{//produce data frame report for R reporting

	LogFile::getInstance()->writeMessage("Starting C Aplpha distance map");
	
	stringstream report;
	report << getHeader() << "\n";
	map<string, Chain*> chains1 = ProteinManager::getInstance()->getChains(pdb1->pdbCode, occupant);
	map<string, Chain*> chains2 = ProteinManager::getInstance()->getChains(pdb2->pdbCode, occupant);
	for (map<string, Chain*>::iterator iter = chains1.begin(); iter != chains1.end(); ++iter)
	{
		Chain* chain = iter->second;
		string chainid = iter->first;
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb1->pdbCode, occupant, chainid);
		vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb1->pdbCode, occupant, chainid);

		for (unsigned int i = 0; i < calphas.size(); ++i)
		{
			Atom* a = calphas[i];
			AminoAcid* aa = aminos[a->aminoId];
			stringstream ss;
			ss << "Distances for: chain " << chain->chainId << ":" << a->atomId << ":" << a->aminoId << ":" << aa->aminoCode;
			LogFile::getInstance()->writeMessage(ss.str());

			for (map<string, Chain*>::iterator iterb = chains2.begin(); iterb != chains2.end(); ++iterb)
			{
				Chain* chainb = iterb->second;
				vector<Atom*> atomsb = chainb->getCAlphas();
				map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb2->pdbCode, occupant, chainb->chainId);
				unsigned int start = 0;								
				for (unsigned int j = start; j < atomsb.size(); ++j)
				{
					Atom* b = atomsb[j];
					AminoAcid* ab = aminosb[b->aminoId];
					double distance = a->atomicDistance(b,shifted);
					if (distance < 25) //TODO this should be configurable					
						report << getRow(aa, a, ab, b, distance,false) << "\n";											
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

string CAlphaReport::getRow(AminoAcid * aa, Atom* a, AminoAcid* ab, Atom* b, double distance, bool singleChain)
{
	int aid = aa->structureAaId;
	int bid = ab->structureAaId;
	if (singleChain)
	{
		aid = aa->aminoId;
		bid = ab->aminoId;
	}
	
	//The amino ID needs to be its place in the entire structure.
	stringstream report;
	report << a->aminoCode << "," << aid << ",";
	report << aa->chainId << ",";
	report << aa->getSS() << ",";
	report << aa->Hydro << ",";
	report << aa->Donicity << ",";
	report << aa->Chemical << ",";
	report << aa->Polar << ",";
	report << b->aminoCode << "," << bid << ",";
	report << ab->chainId << ",";
	report << ab->getSS() << ",";
	report << ab->Hydro << ",";
	report << ab->Donicity << ",";
	report << ab->Chemical << ",";
	report << ab->Polar << ",";
	report << distance;
	return report.str();
}

string CAlphaReport::getHeader()
{
	return "Amino1,Id1,Chain1,SS1,Hydro1,Donicity1,Chemical1,Polar1,Amino2,Id2,Chain2,SS2,Hydro2,Donicity2,Chemical2,Polar2,Distance";
}

