#include "AminoAcid.h"
#include "ProteinManager.h"

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

void AminoAcid::createBonds(AminoAcid* aaP, AminoAcid* aaPP)
{ //Bonds are formed in the order
	//C'-N-CA-C-N''-CA''
	_aaPrev = aaP;
	_aaNext = aaPP;

	//Set all atom references
	_atmCp = _aaPrev->_atoms["C"];
	_atmCApp = _aaNext->_atoms["CA"];;
	_atmNpp = _aaNext->_atoms["N"];;

	_torsion = new BackboneTorsion(aminoCode, aminoId, _atmCp, _atoms["N"], _atoms["CA"], _atoms["C"], _atmNpp, _atmCApp);

	//	_sideTorsion = createSideChainTorsion();
	//TODO Set bonds too but this could be done by looking at the chain info automatically

}

void AminoAcid::add(Atom* atom)
{
	_atoms.insert(pair<string, Atom*>(atom->elementName, atom));
}
