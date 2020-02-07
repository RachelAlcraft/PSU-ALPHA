#include "Chain.h"

Chain::Chain(string pdb_code,string chain_id)
{
	pdbCode = pdb_code;
	chainId = chain_id;
	dataId = pdbCode;
}

Chain::~Chain() //responsible for amino acids
{
	for (map<int, AminoAcid*>::iterator iter = _aminos.begin(); iter != _aminos.end(); ++iter)
	{
		delete iter->second;
	}
	_aminos.clear();
}

AminoAcid* Chain::getAminoAcid(int aminoId)
{
	map<int,AminoAcid*>::iterator iter = _aminos.find(aminoId);
	if (iter == _aminos.end())	
		return nullptr;
	else	
		return iter->second;		
}

vector<Atom*> Chain::getCAlphas()
{
	vector<Atom*> calphas;
	for (map<int, AminoAcid*>::iterator iter = _aminos.begin(); iter != _aminos.end(); ++iter)
	{
		Atom* calpha = iter->second->getCAlpha();
		if (calpha)
			calphas.push_back(calpha);
	}
	return calphas;
}

void Chain::addAminoAcid(AminoAcid* aa)
{
	map<int, AminoAcid*>::iterator iter = _aminos.find(aa->aminoId);
	if (iter == _aminos.end())
		_aminos.insert(pair<int, AminoAcid*>(aa->aminoId,aa));

}
