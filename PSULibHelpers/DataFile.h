#pragma once
#include <string>
#include <vector>
/*
This prints out a file in a data format for R.
Comma delim with headers
Each row should be an observation
*/

using namespace std;

class DataFile
{
private:
	string _filename;
public:
	vector<string> headerVector;
	vector<vector<string>> fileVector;
public:
	DataFile(string filepath);
	void print();
};

