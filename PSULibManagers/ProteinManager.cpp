#include "ProteinManager.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include "PDBFile.h"
#include <LogFile.h>
#include <StringManip.h>


using namespace std;

ProteinManager* ProteinManager::instance = 0;

ProteinManager::ProteinManager()
{

}
ProteinManager::~ProteinManager() //responsible for pdbs
{
	for (map<string, PDBFile*>::iterator iter = _pdbfiles.begin(); iter != _pdbfiles.end(); ++iter)
	{
		delete iter->second;
	}
	_pdbfiles.clear();
}

ProteinManager* ProteinManager::getInstance()
{
	if (!instance)
		instance = new ProteinManager();
	return instance;
}

void ProteinManager::createConfigData(string path)
{
	_configPath = path;
	//AMINO ACID GENERAL DATA//////////////////////////////////////////////
	//Load the amino acid general info, note the chi definitions are actually used to create chi angles
	string aafile = path + "data_aminoinfo.csv";
	vector<string> file;
	ifstream myfile(aafile);
	if (myfile.is_open())
	{
		string line = "";
		while (getline(myfile, line))
			file.push_back(line);
	}
	//Each vector starts with 3 letter code and then has a string of info, comma delim
	for (unsigned int i = 0; i < file.size(); ++i)
	{
		string line = file[i];
		int pos = line.find(",");
		if (pos > 0)
		{
			string aa = line.substr(0, pos);			
			if (aa.length() == 3)
				_aa_generaldata.insert(pair<string, string>(aa, line));//no error checking! It is fixed data. Should be correct!
		}
	}


	//GEOMETRIC VALUES//////////////////////////////////////////////

	//FORCE FIELD///////////////////////////////////////////////////

}

PDBFile* ProteinManager::getOrAddPDBFile(string pdbCode, string filename)
{
	map<string, PDBFile*>::iterator iter = _pdbfiles.find(pdbCode);
	if (iter == _pdbfiles.end())
	{
		PDBFile* pdb = new PDBFile(filename,pdbCode);
		_pdbfiles.insert(pair<string, PDBFile*>(pdbCode, pdb));
		return pdb;
	}
	else
	{
		return iter->second;
	}
}

vector<string> ProteinManager::getAminoData(string aminoCode)
{
	return StringManip::stringToVector(_aa_generaldata[aminoCode],",");//unvalidated. Whole system relies on data being good for now.
}

void ProteinManager::addAtom(string pdbCode, string chainId, int aminoId, Atom* atm)
{
	PDBFile* pdb = _pdbfiles[pdbCode];
	Chain* chain = pdb->getChain(chainId);
	AminoAcid* aa = chain->getAminoAcid(aminoId);
	aa->add(atm);
}


AminoAcid* ProteinManager::getOrAddAminoAcid(string pdbCode, string chainId, int aminoId, string aminoCode)
{
	stringstream id;
	id << pdbCode << chainId << aminoId;

	PDBFile* pdb = _pdbfiles[pdbCode];
	Chain* chain = pdb->getChain(chainId);
	if (chain)
	{
		AminoAcid* aa = chain->getAminoAcid(aminoId);		
		if (!aa)
		{
			AminoAcid* aa = new AminoAcid(pdbCode, chainId, aminoId, aminoCode);
			chain->addAminoAcid(aa);
			return aa;
		}
		else
		{
			return aa;
		}
	}
	else
	{
		LogFile::getInstance()->writeMessage("Amino acid not found " + pdbCode);
		return nullptr;
	}
	
}
Chain* ProteinManager::getOrAddChain(string pdbCode, string chainId)
{
	PDBFile* pdb = _pdbfiles[pdbCode];
	Chain* chain = pdb->getChain(chainId);
	if (!chain)
	{
		Chain* ch = new Chain(pdbCode, chainId);
		pdb->addChain(ch);
		return ch;
	}
	else
	{
		return chain;
	}	
}
map<string, Chain*> ProteinManager::getChains(string pdbCode)
{
	PDBFile* pdb = _pdbfiles[pdbCode];
	return pdb->getChains();
}
map<int, AminoAcid*> ProteinManager::getAminoAcids(string pdbCode, string chainId)
{
	PDBFile* pdb = _pdbfiles[pdbCode];
	Chain* chain = pdb->getChain(chainId);
	return chain->getAminoAcids();

}

vector<Atom*>  ProteinManager::getCAlphas(string pdbCode, string chainId)
{
	PDBFile* pdb = _pdbfiles[pdbCode];
	Chain* chain = pdb->getChain(chainId);
	return chain->getCAlphas();	
}




