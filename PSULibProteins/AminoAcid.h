#pragma once
/*
It is annoying to work with 20 inherited subclasses
So this is an effort to control all the indivudual aa data from a data file
This includes the definitons of the torsions angles
Each aa will have its general information and its specific information, which is basically its atoms
*/
#include <string>
#include <map>
#include <Torsion.h>
#include "Atom.h"

using namespace std;


class AminoAcid
{
public:
	//General information
	string AminoCode;
	string AminoName;	
	string AminoLetter;
	double Hydro;
	string Hydropathy;
	double Volume;
	string Donicity;
	string Chemical;
	string Physio;
	int Charge;
	bool Polar;
	string Formula;
	string AtomChain;
	string Chi1;
	string Chi2;
	string Chi3;
	string Chi4;
	string Chi5;

	//Specific information
	//Unique
	string dataId;
	int aminoId;//unique id within the pdb file
	int structureAaId;//this is within the structure within the chain (for contact maps)
	string aminoCode; //the 3 letter code
	AminoAcid* _aaPrev;
	AminoAcid* _aaNext;
	//For chi,psi,omega
	Atom* _atmCp;
	Atom* _atmNpp;
	Atom* _atmCApp;	
	BackboneTorsion* _torsion;
	SidechainTorsion* _sideTorsion;
	map<string, Atom*> _atoms;
	//Parent
	string pdbCode;
	string chainId;

	//Child
	//do the atoms need to be here or ok just in the manager?
	vector<AtomBond> _bonds;
	vector<AtomAngle> _angles;
	vector<AtomTorsion> _torsions;

public:
	AminoAcid(string pdb_code, string chain_id, int amino_id, int structure_id, string aminoCode);
	~AminoAcid(); //responsible for the torsions
	void createBonds(AminoAcid* aaP, AminoAcid* aaPP);
	void add(Atom*);
	BackboneTorsion* getBackboneTorsion() { return _torsion; }
	SidechainTorsion* getSidechainTorsion() { return _sideTorsion; }
	Atom* getCAlpha() { return _atoms["CA"]; }//obviously this needs error checking, but the whole system fails if there are no CA on each aa so it doesn't matter how badly it crashes TODO
	map<string, Atom*> getAtoms() { return _atoms; }
	string getSS();
	void createScoringAtoms();
	vector<AtomBond> getAtomBonds() { return _bonds; }
	vector<AtomAngle> getAtomAngles() { return _angles; }
	vector<AtomTorsion> getAtomTorsions() { return _torsions; }
//New GeoDefinitions
	vector<AtomGeo*> getAtomBonds(vector<string> atoms);
	vector<AtomGeo*> getAtomCAlphas(vector<string> atoms);
	vector<AtomGeo*> getAtomOneFours(vector<string> atoms);
	vector<AtomGeo*> getAtomInter(vector<string> atoms);
	vector<AtomGeo*> getAtomAngles(vector<string> atoms);
	vector<AtomGeo*> getAtomDihedrals(vector<string> atoms);
	vector<AtomGeo*> getAtomImpropers(vector<string> atoms);

private: //Helper functions
	vector<Atom*> atomsFromString(string atomstring); //we pass around atom strings like N-CA-C and would like to be able to turn those into acrual atoms
	

};


