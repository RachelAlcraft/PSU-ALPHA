#pragma once

#include <string>
#include <vector>
#include <map>
#include "AminoAcid.h"
#include <GeoTransformation.h>
//#include "Bond.h"
#include "Atom.h"
#include "Chain.h"

class PDBFile
{
private:
	string _filename;
	vector<string> _file; //The original file is saved so we can print it back out
public:
	string pdbCode;
	

private:
	map<string, Chain*> _chains;

public:
	PDBFile(string,string);
	~PDBFile();
	string getFileString();	
	void addLinks();
	map<string, Chain*> getChains() { return _chains; }
	Chain* getChain(string chainId);
	void addChain(Chain* ch);
	void loadData();
	void applyTransformation(GeoTransformation trans);


private:
	void createFileVector();
	
	

};

