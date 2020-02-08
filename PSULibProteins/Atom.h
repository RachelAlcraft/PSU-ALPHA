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
	double atomicDistance(Atom*);
	void applyTransformation(GeoTransformation trans);

private:
	//string trim(string);

};

