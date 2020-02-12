#include "CAlphaReport.h"
#include <LogFile.h>
#include <ProteinManager.h>

void CAlphaReport::printReport(PDBFile* pdb, string chain1, string chain2, string fileName)
{//produce data frame report for R reporting
	
	LogFile::getInstance()->writeMessage("Starting C Aplpha distance map");
	if (chain1 != "" && chain2 != "")
	{
		printSingleChainReport(pdb, chain1, chain2, fileName);
	}
	else
	{
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
						double distance = a->atomicDistance(b,false);
						if (distance < 25) //TODO this should be configurable
						{
							report << getRow(aa, a, ab, b, distance) << "\n";							
							if (!(chainself & (i == j)))//and reverse to halve the time														
								report << getRow(ab, b, aa, a, distance) << "\n";							
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

void CAlphaReport::printSingleChainReport(PDBFile* pdb, string chain1, string chain2, string fileName)
{//produce data frame report for R reporting

	stringstream report;
	report << "Amino1, Chain1, Hydro1, Donicity1, Chemical1, Polar1, Amino2, Chain2, Hydro2, Donicity2, Chemical2, Polar2, Distance\n";
	Chain* pchain1 = ProteinManager::getInstance()->getOrAddChain(pdb->pdbCode,chain1);
	Chain* pchain2 = ProteinManager::getInstance()->getOrAddChain(pdb->pdbCode, chain2);		
	map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, chain1);
	vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb->pdbCode, chain1);

	for (unsigned int i = 0; i < calphas.size(); ++i)
	{
		Atom* a = calphas[i];
		AminoAcid* aa = aminos[a->aminoId];
		stringstream ss;
		ss << "Distances for: chain " << chain1 << ":" << a->atomId << ":" << a->aminoId << ":" << aa->aminoCode;
		LogFile::getInstance()->writeMessage(ss.str());
				
		vector<Atom*> atomsb = pchain2->getCAlphas();
		map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb->pdbCode, chain2);				
		for (unsigned int j = 0; j < atomsb.size(); ++j)
		{
			Atom* b = atomsb[j];
			AminoAcid* ab = aminosb[b->aminoId];
			double distance = a->atomicDistance(b,false);
			if (distance < 25) //TODO this should be configurable
			{				
				report << getRow(aa, a, ab, b, distance) << "\n";
			}
		}			
	}
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		outfile << report.str();
	}
	
}

void CAlphaReport::printMultiReport(PDBFile* pdb1, PDBFile* pdb2, string fileName, bool shifted)
{//produce data frame report for R reporting

	LogFile::getInstance()->writeMessage("Starting C Aplpha distance map");
	
	stringstream report;
	report << "Amino1,PDB1,Hydro1,Donicity1,Chemical1,Polar1,Amino2,PDB2,Hydro2,Donicity2,Chemical2,Polar2,Distance\n";
	map<string, Chain*> chains1 = ProteinManager::getInstance()->getChains(pdb1->pdbCode);
	map<string, Chain*> chains2 = ProteinManager::getInstance()->getChains(pdb2->pdbCode);
	for (map<string, Chain*>::iterator iter = chains1.begin(); iter != chains1.end(); ++iter)
	{
		Chain* chain = iter->second;
		string chainid = iter->first;
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdb1->pdbCode, chainid);
		vector<Atom*> calphas = ProteinManager::getInstance()->getCAlphas(pdb1->pdbCode, chainid);

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
				map<int, AminoAcid*> aminosb = ProteinManager::getInstance()->getAminoAcids(pdb2->pdbCode, chainb->chainId);
				unsigned int start = 0;								
				for (unsigned int j = start; j < atomsb.size(); ++j)
				{
					Atom* b = atomsb[j];
					AminoAcid* ab = aminosb[b->aminoId];
					double distance = a->atomicDistance(b,shifted);
					if (distance < 25) //TODO this should be configurable					
						report << getRow(aa, a, ab, b, distance) << "\n";											
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

string CAlphaReport::getRow(AminoAcid * aa, Atom* a, AminoAcid* ab, Atom* b, double distance)
{
	//The amino ID needs to be its place in the entire structure.
	stringstream report;
	report << a->aminoCode << "," << aa->structureAaId << ",";
	report << aa->Hydro << ",";
	report << aa->Donicity << ",";
	report << aa->Chemical << ",";
	report << aa->Polar << ",";
	report << b->aminoCode << "," << ab->structureAaId << ",";
	report << ab->Hydro << ",";
	report << ab->Donicity << ",";
	report << ab->Chemical << ",";
	report << ab->Polar << ",";
	report << distance;
	return report.str();
}

