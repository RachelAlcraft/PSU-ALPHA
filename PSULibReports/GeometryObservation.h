#pragma once
#include <string>
#include <vector>

using namespace std;


class GeometryObservation
{
	//PdbCode,Chain,AminoAcid,AminoNo,PdbAtoms,SecStruct,GeoType,ExperimentalMethod,GeoAtoms,Value
public:
	double value;
	string pdbCode;
	string aminoCode;
	string aminoNo;
	string pdbAtoms;
	string secStruc;
	string geoType;
	string geoAtoms;
	
};

