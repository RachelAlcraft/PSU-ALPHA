#include "PDBFile.h"
#include <fstream>
#include <iostream>
#include "Atom.h"
#include "ProteinManager.h"
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "LogFile.h"

using namespace std;

PDBFile::PDBFile(string filename, string pdb_code)
{
	pdbCode = pdb_code;	
	_filename = filename;
}

PDBFile::~PDBFile()
{
}

void PDBFile::loadData()
{	
	createFileVector();
	Chain* currentChain = NULL;
	//AminoAcid* currentAmino = NULL;	
	for (unsigned int i = 0; i < _file.size(); ++i)
	{
		string line = _file[i];
		int pos = line.find("ATOM");
		if (pos == 0)
		{
			Atom* atom = new Atom(pdbCode,line);
			if (i % 500 == 0)
				atom->printAtom();
			int id = atom->atomId;
			
			string amino = atom->aminoCode;
			int amino_id = atom->aminoId;
			string chain = atom->chainId;
			Chain* pchain = ProteinManager::getInstance()->getOrAddChain(pdbCode, chain);
			AminoAcid* paa = ProteinManager::getInstance()->getOrAddAminoAcid(pdbCode, chain, amino_id, amino);
			ProteinManager::getInstance()->addAtom(pdbCode,chain,amino_id,atom);															
		}
	}
	addLinks();
}

string PDBFile::getFileString()
{
	string wholefile = "";
	ifstream myfile(_filename);
	if (myfile.is_open())
	{
		string line = "";
		while (getline(myfile, line))
			wholefile += line;
	}
	return wholefile;
}

void PDBFile::createFileVector()
{

	ifstream myfile(_filename);
	if (myfile.is_open())
	{
		string line = "";
		while (getline(myfile, line))
			_file.push_back(line);
	}
}

void PDBFile::addLinks()
{
	map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdbCode);
	for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
	{
		Chain* chain = iter->second;		
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdbCode, chain->chainId);
		AminoAcid* lastaa = NULL;
		AminoAcid* nextaa = NULL;
		map<int, AminoAcid*>::iterator iteraa = aminos.begin();
		while (iteraa != aminos.end())
		{
			AminoAcid* aa = iteraa->second;
			++iteraa;
			if (iteraa != aminos.end())
			{
				nextaa = iteraa->second;
				if (lastaa && iteraa != aminos.end())//For now ignore first and last of each chain as per chimera ramachandran plots
				{
					aa->createBonds(lastaa, nextaa);
				}
			}
			lastaa = aa;
		}
	}
}

Chain* PDBFile::getChain(string chainId)
{
	map<string, Chain*>::iterator iter = _chains.find(chainId);
	if (iter == _chains.end())
		return nullptr;
	else
		return iter->second;
}

void PDBFile::addChain(Chain* ch)
{
	map<string, Chain*>::iterator iter = _chains.find(ch->chainId);
	if (iter == _chains.end())
		_chains.insert(pair<string, Chain*>(ch->chainId, ch));
}

