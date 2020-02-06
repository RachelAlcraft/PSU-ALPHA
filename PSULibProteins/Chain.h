#pragma once

#include <string>
#include <map>
#include <AminoAcid.h>

using namespace std;

class Chain
{
public:
	//Parent
	string pdbCode;
	string dataId;
	string chainId;

	//Children
	map<int, AminoAcid*> _aminos;

public:
	Chain(string pdb_code, string chain_id);
	AminoAcid* getAminoAcid(int aminoId);
	map<int, AminoAcid*> getAminoAcids() { return _aminos; }
	void addAminoAcid(AminoAcid* aa);
};

