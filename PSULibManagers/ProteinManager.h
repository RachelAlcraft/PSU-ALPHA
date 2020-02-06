#pragma once
/* 
SINGLETON pattern datamanager class
This is a kind of runtime relational database
The class also manages the loading of some config data
The important thing about this class it manages all the pointes and is the responsible object for all pointer deletion
Pointers are created and deleted just once in the program not really dynamically
*/
#include <string>
#include <map>
#include <set>
#include <vector>
#include <AminoAcid.h>
#include <Chain.h>
#include <Atom.h>
#include <PDBFile.h>

using namespace std;
class ProteinManager
{
private:
	static ProteinManager* instance;
	ProteinManager();
	string _configPath;

	//DATA HIERARCHY//////////////////////////////////////////
	//Amino acids have general data
	map<string, string> _aa_generaldata; // 3 letter AA code vs a long string of data to parse into the aa init
	
	//Amino acids have geometric preferences

	//Force field

	
	//PROTEIN HIERARCHY////////////////////////////////////
	//a pdb file has
		//chains
			//which have amino acids
				//which have atoms
	//There are also bonds, torsions angles and angles

	//PDBFile has chains
	map<string,PDBFile*> _pdbfiles;
	//map<string, vector<Chain*>> _chains; //id = PDBCODE 		
	//Chains have amino acids
	//map<string, vector<AminoAcid*>> _aminos; //id = PDBCODE_chainNo 
	//amino acids have atoms (which may be shifted in optimisation algorithms)
	//map<string, vector<Atom*>> _atoms; //id = PDBCODE_chainNo_aminoNo

	//amino acids have bonds

	//amino acids have angles
	
	//amino acids have torsion phi, psi, omega and chi angles

	//Amino acids have improper bonds
	
	//Atoms are involved in hydrogen bonds

public:
	static ProteinManager* getInstance();
	void createConfigData(string path);
	vector<string> getAminoData(string aminoCode);
	void addAtom(string pdbCode, string chainId, int aminoId, Atom* atm);
	PDBFile* getOrAddPDBFile(string pdbCode, string filename);
	AminoAcid* getOrAddAminoAcid(string pdbCode, string chainId, int aminoId, string amino_code);		
	Chain* getOrAddChain(string pdbCode, string chainId);
	map<string, Chain*> getChains(string pdbCode);	
	map<int, AminoAcid*> getAminoAcids(string pdbCode, string chainId);		

private:
	vector<string> stringToVector(string input, string delim);






};

