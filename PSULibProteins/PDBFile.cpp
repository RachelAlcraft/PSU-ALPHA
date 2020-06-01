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
#include <NucleicAcid.h>
#include <StringManip.h>

using namespace std;

PDBFile::PDBFile(string filename, string pdb_code)
{
	pdbCode = pdb_code;	
	_filename = filename;
	//residues = 0;
	experimentalMethod = "XR"; // only XR at the moment
	loadedText = false;
	loadedAminos = false;
	loadedBonds = false;
	loadedTorsions = false;
}

PDBFile::~PDBFile()
{
	for (map<string, ProteinStructure*>::iterator iter = _proteinVersions.begin(); iter != _proteinVersions.end(); ++iter)
	{
		delete iter->second;
	}
	_proteinVersions.clear();
}

void PDBFile::loadData()
{	
	createFileVector();	
	structureInfo();
	loadedText = true;
}

void PDBFile::structureInfo()
{
	/*
	This goes through the structure and loads info
	*/
	rvalue = "NULL";
	rfree = "NULL";
	inComplex = "N";
	institution = "";
	nullModel = false;
	
	for (unsigned int i = 0; i < _file.size(); ++i)
	{
		string line = _file[i];
		
		int posP = line.find("3   PROGRAM");		
		int posA = line.find("JRNL        AUTH");		
		int posR = line.find("3   R VALUE");
		int posRF = line.find("3   FREE R VALUE     ");
		int posT = line.find("TITLE");	
		int posM = line.find("NCS MODEL : NULL");
		if (i == 0)//header row has date and class in fixed positions
		{
			proteinclass = StringManip::removeChar(StringManip::trim(line.substr(10, 40)),","," ");
			date = StringManip::trim(line.substr(50, 9));
		}
		else if (posP >= 0)
		{
			software = StringManip::removeChar(StringManip::trim(line.substr(27)),","," ");
		}
		else if (posA >= 0)
		{
			institution += StringManip::trim(line.substr(19));
			institution = StringManip::removeChar(institution, ",","-");			
		}		
		else if (posR >= 0)
		{
			if (rvalue == "NULL")
			{
				vector<string> rvec = StringManip::stringToVector(line, ":");
				if (rvec.size() > 1)
					rvalue = StringManip::trim(rvec[1]);
			}
		}
		else if (posRF >= 0)
		{
			if (rfree == "NULL")
			{
				vector<string> rvec = StringManip::stringToVector(line, ":");
				if (rvec.size() > 1)
					rfree = StringManip::trim(rvec[1]);
			}
		}		
		else if (posT >= 0 && inComplex == "N")
		{
			int pos6 = line.find("COMPLEX");
			inComplex = pos6 > 0 ? "Y" : "N";			
		}
		else if (posM >= 0)
		{
			LogFile::getInstance()->writeMessage("NCS Model=TRUE");
			nullModel = true;
		}
	}	
}


void PDBFile::loadAtoms()
{
	/*
	Clearly this is not how it should be done but
	before we start, not knowing the occupancy and needing everything to be copied
	create A,B,S,D occupant versions and then delete the ones we don;t need	
	*/
	addStructureVersion("A");
	addStructureVersion("B");
	addStructureVersion("C");
	addStructureVersion("D");

	Chain* currentChain = NULL;
	int structureId = 0;
	//AminoAcid* currentAmino = NULL;	
	for (unsigned int i = 0; i < _file.size(); ++i)
	{
		string line = _file[i];
		int pos = line.find("ATOM");
		if (pos == 0)
		{
			Atom* atom = new Atom(pdbCode, line);
			if (i % 500 == 0)
				atom->printAtom();
			int id = atom->atomId;

			string amino = atom->aminoCode;
			int amino_id = atom->aminoId;
			string chain = atom->chainId;
			string ocpnt = atom->occupant;
			vector<string> alloccupants;
			if (ocpnt == "X")
			{
				alloccupants.push_back("A");
				alloccupants.push_back("B");
				alloccupants.push_back("C");
				alloccupants.push_back("D");
			}
			else
			{
				alloccupants.push_back(ocpnt);
			}
			atom->isAmino = true;
			if (ProteinManager::getInstance()->isNucleicAcid(amino))			
				atom->isAmino = false;		

			for (unsigned int op = 0; op < alloccupants.size(); ++op)//the atom needs to be added to all the approrpaite occupant structures
			{
				Chain* pchain = ProteinManager::getInstance()->getOrAddChain(pdbCode, alloccupants[op], chain);

				if (atom->isAmino)
				{
					AminoAcid* paa = ProteinManager::getInstance()->getOrAddAminoAcid(pdbCode, alloccupants[op], chain, amino_id, amino, structureId);
					if (paa != nullptr)
						ProteinManager::getInstance()->addAtom(pdbCode, alloccupants[op], chain, amino_id, atom);
				}
				else
				{
					NucleicAcid* pna = ProteinManager::getInstance()->getOrAddNucleicAcid(pdbCode, alloccupants[op], chain, amino_id, amino, structureId, nucleotides);
				}
			}
		}
	}
	//remove unnecessary occupants, this is obviously not ideal! TODO
	vector<string> occupantList = ProteinManager::getInstance()->occupantList(pdbCode);
	if (std::find(occupantList.begin(), occupantList.end(), "B") == occupantList.end())
		removeStructureVersion("B");
	if (std::find(occupantList.begin(), occupantList.end(), "C") == occupantList.end())
		removeStructureVersion("C");
	if (std::find(occupantList.begin(), occupantList.end(), "D") == occupantList.end())
		removeStructureVersion("D");
	

	loadedAminos = true;
	
}
void PDBFile::loadBonds()
{
	addLinks();
	loadedBonds = true;
	loadedTorsions = true;
	vector<string> seqs = getSequence();
	if (seqs.size() > 1)
	{
		//if the chains are not independently verified and are all the same then keep only the first
		string sequence = "";
		bool differ = false;
		for (unsigned int r = 0; r < seqs.size(); ++r)
		{
			if (r == 0)
				sequence = seqs[r];
			else if (seqs[r] != sequence)
				differ = true;
		}
		if (!differ)
		{
			if (!nullModel)
			{
				LogFile::getInstance()->writeMessage("Identical chains with no NCS setting present");
				removeRepeatedChains();
			}
			else
			{
				LogFile::getInstance()->writeMessage("Identical chains but NCS setting present");
			}
		}
		
	}
}
void PDBFile::loadTorsions()
{//currently done as paert of loadBonds
}



void PDBFile::printShiftedFile(string fileRoot)
{
	string fileName = fileRoot + pdbCode + ".pdb";
	map<int,Atom*> atoms = ProteinManager::getInstance()->getAtomsMap(pdbCode,"A");
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

void PDBFile::prepareStructureVersions()
{
}

void PDBFile::addLinks()
{
	for (map<string, ProteinStructure*>::iterator iter = _proteinVersions.begin(); iter != _proteinVersions.end(); ++iter)
	{
		string occupant = iter->first;
		ProteinStructure* ps = iter->second;		
		map<string, Chain*> chains = ps->getChains();
		for (map<string, Chain*>::iterator iter = chains.begin(); iter != chains.end(); ++iter)
		{
			Chain* chain = iter->second;
			map<int, AminoAcid*> aminos = ProteinManager::getInstance()->getAminoAcids(pdbCode, occupant,chain->chainId);
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
						//aa->createBonds(lastaa, nextaa);
						aa->createBonds(lastaa, nextaa);
						aa->createScoringAtoms();
					}
				}

				lastaa = aa;

			}
		}
	}
}

void PDBFile::removeRepeatedChains()
{
	for (map<string, ProteinStructure*>::iterator iter = _proteinVersions.begin(); iter != _proteinVersions.end(); ++iter)
	{
		string occupant = iter->first;
		ProteinStructure* ps = iter->second;
		ps->removeRepeatedChains();		
	}
}

ProteinStructure* PDBFile::getStructureVersion(string occupant)
{

	map<string, ProteinStructure*>::iterator iter = _proteinVersions.find(occupant);
	if (iter == _proteinVersions.end())
		return nullptr;
	else
		return iter->second;
}

void PDBFile::addStructureVersion(string occupant)
{	
	map<string, ProteinStructure*>::iterator iter = _proteinVersions.find(occupant);
	if (iter == _proteinVersions.end())
	{		
		_proteinVersions.insert(pair<string, ProteinStructure*>(occupant, new ProteinStructure(pdbCode, occupant)));		
	}
}

void PDBFile::removeStructureVersion(string occupant)
{
	map<string, ProteinStructure*>::iterator iter = _proteinVersions.find(occupant);
	if (iter != _proteinVersions.end())
	{
		_proteinVersions.erase(occupant);		
	}
}

vector<string> PDBFile::getSequence()
{
	vector<string> sequences;
	ProteinStructure* ps = getStructureVersion("A");
	map<string, Chain*> chains = ps->getChains();
	map<string, Chain*>::iterator iterch = chains.begin();
	for (iterch; iterch != chains.end(); ++iterch)
	{
		string chainseq = "";
		map<int, AminoAcid*> aminos = iterch->second->getAminoAcids();
		map<int, AminoAcid*>::iterator iteraa = aminos.begin();
		for (iteraa; iteraa != aminos.end(); ++iteraa)
		{
			chainseq += iteraa->second->AminoLetter;
		}
		sequences.push_back(chainseq);
	}
	return sequences;
}

string PDBFile::maxChain()
{	
	string letter = "";
	ProteinStructure* ps = getStructureVersion("A");
	map<string, Chain*> chains = ps->getChains();
	map<string, Chain*>::iterator iterch = chains.end();
	--iterch;
	if (iterch != chains.end())
		letter = iterch->first;	
	return letter;
}


