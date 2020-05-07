#pragma once
#include <PDBFile.h>
#include <vector>

using namespace std;

class GeometricalDataReport
{
private:
	/*
	the file is of the format
	-------------------------------
	AminoAcid, GeoType, Atoms
	#Previousand next atoms, ,
	*,BOND,CP-N
	ILE,ONEFOUR,CD1-CG2
	-------------------------------
	Where a # means ignore, and a * means all amino acids
	*/
	vector < vector<string>> geoDefinitions;
	
public:
	GeometricalDataReport(vector < vector<string>> fileVector) { geoDefinitions = fileVector; }
	void printReport(PDBFile* pdb, string filename1/*, string filename2*/, string directory, string geodir);
private:
	void printOneReport(PDBFile* pdb, string occupant, string filename1);// , string filename2);
	void printOneReportWithGeoDef(PDBFile* pdb, string occupant, string fileName1);
	vector<string> getGeoDefinitions(string aminoCode, string geoType);
	string getReportString(AtomGeo* ab, PDBFile* pdb, string occupant, string geoType);
	double quickRound(double val);
};

