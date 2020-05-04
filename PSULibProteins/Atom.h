#pragma once
#include "GeoCoords.h"
#include "GeoVector.h"
#include <string>
#include <vector>
#include <GeoTransformation.h>


using namespace std;
class Atom
{
public: //public struct
	//unique
	string dataId;
	int atomId;
	string elementName;
	string elementType;
	string occupant;
	double occupancy;
	double bfactor;
	bool isAmino;
	
	//Parent
	string pdbCode;	
	string chainId;
	string aminoCode;
	int aminoId;
	
	//Geometric Info	
	GeoCoords coords;
	GeoCoords shifted_coords; // alternatively I could store a vector of shifts which could be animated - TODO
private:
	vector<Atom*> _bonds;
public:
	Atom(string,string);
	~Atom();
	void applyShift(double, double, double, bool);
	void printAtom();
	string getDescription();
	GeoVector vectorDifference(Atom*);
	double atomicDistance(Atom*,bool shifted);
	void applyTransformation(GeoTransformations* trans);

private:
	//string trim(string);

};

class AtomGeo
{
protected:
	Atom* _A1;
	Atom* _A2;
	string _atomString;
	string _atomNoString;
	string _SS;
public:
	string geoDef;
	string allAAs;
public:
	AtomGeo(string ss, Atom* a1, Atom* a2, string geo);
	string getChain() { return _A1->chainId; }
	string getAA() { return _A1->aminoCode; }
	int getId() { return _A1->atomId; }
	string getAtoms() { return _atomString; }
	string getAtomNos() { return _atomNoString; }
	string getSS() { return _SS; }
	virtual double getValue() = 0;
};

class AtomBond: public AtomGeo
{
public:	
	AtomBond(string ss, Atom* a1, Atom* a2, string geo);
	double getValue() override;
};

class AtomAngle : public AtomGeo
{	
	Atom* _A3;
public:
	AtomAngle(string ss, Atom* a1, Atom* a2, Atom* a3, string geo);
	double getValue() override;
	

};

class AtomTorsion: public AtomGeo
{	
	Atom* _A3;
	Atom* _A4;
public:
	AtomTorsion(string ss, Atom* a1, Atom* a2, Atom* a3, Atom* a4, string geo);
	double getValue() override;
	
};

