#include "CSVFile.h"
#include <fstream>
#include <StringManip.h>

CSVFile::CSVFile(string filepath, string sep)
{
	_filename = filepath;

	ifstream myfile(_filename);
	if (myfile.is_open())
	{
		string line = "";
		while (getline(myfile, line))
		{
			//decide to delete some data I am not handling TODO!
			vector<string> row = StringManip::stringToVector(line, sep);
			fileVector.push_back(row);
		}
	}	
}
