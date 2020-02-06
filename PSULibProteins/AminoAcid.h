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
	string Chain;
	string Chi1;
	string Chi2;
	string Chi3;
	string Chi4;
	string Chi5;

	//Specific information
	//Unique
	string dataId;
	int aminoId;//unique id within the pdb file
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
public:
	AminoAcid(string pdb_code, string chain_id, int amino_id, string aminoCode);
	void createBonds(AminoAcid* aaP, AminoAcid* aaPP);
	void add(Atom*);
	BackboneTorsion* getBackboneTorsion() { return _torsion; }
	SidechainTorsion* getSidechainTorsion() { return _sideTorsion; }

};


