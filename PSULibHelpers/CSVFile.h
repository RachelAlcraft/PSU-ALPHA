#pragma once

#include <string>
#include <vector>

using namespace std;

class CSVFile
{
	private: 
	string _filename;
public:
	vector < vector<string>> fileVector;
public:
	CSVFile(string filepath, string sep);
};

