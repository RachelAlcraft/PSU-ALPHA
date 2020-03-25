#include "DataFile.h"
#include <sstream>
#include <fstream>

using namespace std;

DataFile::DataFile(string filepath)
{
	_filename = filepath;
}

void DataFile::print()
{
	stringstream ss;
	for (unsigned int i = 0; i < headerVector.size(); ++i)
	{
		if (i > 0)
			ss << ",";
		ss << headerVector[i];
	}

	for (unsigned int i = 0; i < fileVector.size(); ++i)
	{
		vector<string> row = fileVector[i];
		ss << "\n";
		for (unsigned int j = 0; j < row.size(); ++j)
		{
			if (j > 0)
				ss << ",";
			ss << row[j];
		}		
	}


	ofstream outfile(_filename);
	if (outfile.is_open())
	{
		outfile << ss.str();
	}
}
