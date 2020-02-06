#pragma once
#include "Atom.h"

class Torsion
{//FOR NOW insist that a torsion has all of the possible atoms
	//Bonds are formed in the order
	//C'-N-CA-C-N''-CA''
public:
	int id;
	string aa;
protected:
	Atom* _N;
	Atom* _CA;
	Atom* _C;

public: //methods
	Torsion(string, int, Atom*, Atom*, Atom*);
	double getDihedralAngle(Atom*, Atom*, Atom*, Atom*);

};

class BackboneTorsion : public Torsion
{
protected:
	Atom* _Cp;
	Atom* _Npp;
	Atom* _CApp;
public:
	BackboneTorsion(string, int,Atom*,Atom*,Atom*,Atom*,Atom*,Atom*);
	double getPhi();
	double getPsi();
	double getOmega();

};

class SidechainTorsion : public Torsion
{
protected:
	Atom* _AG1;
	Atom* _AD1;

public:
	SidechainTorsion(string, int, Atom*, Atom*, Atom*, Atom*, Atom*);
	double getChi1();
	double getChi2();

};

