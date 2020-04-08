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
	double occupancy;
	
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

class AtomBond
{
protected:
	Atom* _A1;
	Atom* _A2;
	string _atomString;
	string _SS;
public:
	AtomBond(string ss, Atom* a1, Atom* a2);
	string getChain() { return _A1->chainId; }
	string getAA() { return _A1->aminoCode; }
	int getId() { return _A1->atomId; }	
	string getAtoms() { return _atomString; }
	string getSS() { return _SS; }
	virtual double getValue();
};

class AtomAngle : public AtomBond
{	
	Atom* _A3;
public:
	AtomAngle(string ss, Atom* a1, Atom* a2, Atom* a3);
	double getValue() override;

};

class AtomTorsion: public AtomBond
{	
	Atom* _A3;
	Atom* _A4;
public:
	AtomTorsion(string ss, Atom* a1, Atom* a2, Atom* a3, Atom* a4);
	double getValue() override;
};

