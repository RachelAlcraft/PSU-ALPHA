#pragma once
/*
It is annoying to work with 20 inherited subclasses
So this is an effort to control all the indivudual aa data from a data file
This includes the definitons of the torsions angles
Each aa will have its general information and its specific information, which is basically its atoms
*/
#include <string>

using namespace std;


class AminoAcid
{
public:
	//General information
	string AminoName;
	string AminoCode;
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

	//Parent
	string pdb_code;
	int chain;

	//Child
	//do the atoms need to be here or ok just in the manager?
};

