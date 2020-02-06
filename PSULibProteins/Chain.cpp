#include "Chain.h"

Chain::Chain(string pdb_code,string chain_id)
{
	pdbCode = pdb_code;
	chainId = chain_id;
	dataId = pdbCode;
}

AminoAcid* Chain::getAminoAcid(int aminoId)
{
	map<int,AminoAcid*>::iterator iter = _aminos.find(aminoId);
	if (iter == _aminos.end())	
		return nullptr;
	else	
		return iter->second;		
}

void Chain::addAminoAcid(AminoAcid* aa)
{
	map<int, AminoAcid*>::iterator iter = _aminos.find(aa->aminoId);
	if (iter == _aminos.end())
		_aminos.insert(pair<int, AminoAcid*>(aa->aminoId,aa));

}
