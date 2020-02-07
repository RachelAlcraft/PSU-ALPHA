#include "AminoAcid.h"
#include "ProteinManager.h"
#include <StringManip.h>

AminoAcid::AminoAcid(string pdb_code, string chain_id, int amino_id, string amino_Code)
{
	aminoId = amino_id;	
	pdbCode = pdb_code;
	chainId = chain_id;
	aminoCode = amino_Code;
	vector<string> aa_general = ProteinManager::getInstance()->getAminoData(aminoCode);
	//general data setttings
	AminoCode = aa_general[0];
	AminoName = aa_general[1];
	AminoLetter = aa_general[2];
	Hydro = atof(aa_general[3].c_str());
	Hydropathy = aa_general[4];
	Volume = atof(aa_general[5].c_str());
	Donicity = aa_general[6];
	Chemical = aa_general[7];
	Physio = aa_general[8];
	Charge = atol(aa_general[9].c_str());
	Polar = aa_general[10]=="T"?true:false;
	Formula = aa_general[11];
	Chain = aa_general[12];
	Chi1 = aa_general[13];
	Chi2 = aa_general[14];
	Chi3 = aa_general[15];
	Chi4 = aa_general[16];
	Chi5 = aa_general[17];
}

AminoAcid::~AminoAcid()
{
	delete _torsion;
	delete _sideTorsion;
	_torsion = nullptr;
	_sideTorsion = nullptr;
	_aaPrev = nullptr;
	_aaNext = nullptr;
	_atmCp = nullptr;
	_atmNpp = nullptr;
	_atmCApp = nullptr;
	
	for (map<string, Atom*>::iterator iter = _atoms.begin(); iter != _atoms.end(); ++iter)
	{
		delete iter->second;
	}
	_atoms.clear();
	
}

void AminoAcid::createBonds(AminoAcid* aaP, AminoAcid* aaPP)
{ //Bonds are formed in the order
	//C'-N-CA-C-N''-CA''
	_aaPrev = aaP;
	_aaNext = aaPP;

	//Set all atom references
	_atmCp = _aaPrev->_atoms["C"];
	_atmCApp = _aaNext->_atoms["CA"];
	_atmNpp = _aaNext->_atoms["N"];

	vector<Atom*> backbones = atomsFromString("CP-N-CA-C-NPP-CAPP");
	_torsion = new BackboneTorsion(aminoCode, aminoId, backbones);

	vector<vector<Atom*>> sidechains;
	if (Chi1 != "None")
		sidechains.push_back(atomsFromString(Chi1));
	if (Chi2 != "None")
		sidechains.push_back(atomsFromString(Chi2));
	if (Chi3 != "None")
		sidechains.push_back(atomsFromString(Chi3));
	if (Chi4 != "None")
		sidechains.push_back(atomsFromString(Chi4));
	if (Chi5 != "None")
		sidechains.push_back(atomsFromString(Chi5));
	
	_sideTorsion = new SidechainTorsion(aminoCode, aminoId, sidechains);

	//TODO Set bonds too but this could be done by looking at the chain info automatically
	//AND angles which we will need at some point

}

void AminoAcid::add(Atom* atom)
{
	_atoms.insert(pair<string, Atom*>(atom->elementName, atom));
}

vector<Atom*> AminoAcid::atomsFromString(string atomstring)
{//For beginning and end of chains we return nothing by checking the pointers to the previous and next
	vector<string> stratoms = StringManip::stringToVector(atomstring, "-");
	vector<Atom*> vecatoms;
	bool bCompleteSet = true;
	for (unsigned int i = 0; i < stratoms.size(); ++i)
	{
		Atom* atm = nullptr;
		string stratom = stratoms[i];
		if (stratom == "CP")
		{
			if (_atmCp)
				atm = _atmCp;
			else
				bCompleteSet = false;
		}
		else if (stratom == "NPP")
		{
			if (_atmNpp)
				atm = _atmNpp;
			else
				bCompleteSet = false;
		}
		else if (stratom == "CAPP")
		{
			if (_atmCApp)
				atm = _atmCApp;
			else
				bCompleteSet = false;
		}
		else
			atm = _atoms[stratom];
		if (atm)
			vecatoms.push_back(atm);
	}
	if (bCompleteSet)
		return vecatoms;
	else
		return vector<Atom*>();
}
