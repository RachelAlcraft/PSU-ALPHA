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
	int structureId = 0;
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
			AminoAcid* paa = ProteinManager::getInstance()->getOrAddAminoAcid(pdbCode, chain, amino_id, amino, structureId);
			ProteinManager::getInstance()->addAtom(pdbCode,chain,amino_id,atom);															
		}
	}
	addLinks();
}

void PDBFile::applyTransformation(GeoTransformations* trans)
{
	vector<Atom*> atoms = ProteinManager::getInstance()->getAtoms(pdbCode);
	for (unsigned int i = 0; i < atoms.size(); ++i)
	{
		if (atoms[i])
			atoms[i]->applyTransformation(trans);
	}
}

void PDBFile::printShiftedFile(string fileRoot)
{
	string fileName = fileRoot + pdbCode + ".pdb";
	map<int,Atom*> atoms = ProteinManager::getInstance()->getAtomsMap(pdbCode);
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		for (unsigned int i = 0; i < _file.size(); ++i)
		{
			string line = _file[i];
			int pos = line.find("ATOM");
			if (pos == 0)
			{
				Atom* dummy = new Atom(pdbCode,line);				
				int id = dummy->atomId;
				map<int, Atom*>::iterator atom_iter = atoms.find(id);
				if (atom_iter != atoms.end())
				{
					//31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
					//39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
					//47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
					stringstream ssx, ssy, ssz;
					ssx << setprecision(3) << fixed << atom_iter->second->shifted_coords.x;
					ssy << setprecision(3) << fixed << atom_iter->second->shifted_coords.y;
					ssz << setprecision(3) << fixed << atom_iter->second->shifted_coords.z;
					while (ssx.str().length() < 8)
						ssx << " ";
					while (ssy.str().length() < 8)
						ssy << " ";
					while (ssz.str().length() < 8)
						ssz << " ";
					line.replace(31, 8, ssx.str());
					line.replace(39, 8, ssy.str());
					line.replace(47, 8, ssz.str());
					outfile << line << '\n';
				}
				else
				{
					outfile << _file[i] << '\n'; //just print out what we have or should I error ? TODO
				}
			}
			else
			{
				outfile << _file[i] << '\n'; //just print out what we have or should I error ? TODO
			}
		}
	}
}

string PDBFile::getSequence() // TODO make it per chain if required
{
	string seq = "";
	map<string, Chain*> chains = ProteinManager::getInstance()->getChains(pdbCode);
	for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
	{
		Chain* chain = iter->second;
		map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdbCode, chain->chainId);
		for (map<int, AminoAcid*>::iterator iteraa = aminos.begin(); iteraa != aminos.end(); ++iteraa)
		{
			seq += iteraa->second->AminoLetter;
		}		
	}
	return seq;
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
		{
			//decide to delete some data I am not handling TODO!
			int pos = line.find("HETATM");
			if (pos < 0)
				_file.push_back(line);

		}
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

map<int, Atom*> PDBFile::getAtoms(string pdbCode)
{
	map<int, Atom*> atoms;
	for (map<string, Chain*>::iterator iter = _chains.begin(); iter != _chains.end(); ++iter)
	{
		Chain* ch = iter->second;

	}
	
	return atoms;
}

void PDBFile::addChain(Chain* ch)
{
	map<string, Chain*>::iterator iter = _chains.find(ch->chainId);
	if (iter == _chains.end())
		_chains.insert(pair<string, Chain*>(ch->chainId, ch));
}

