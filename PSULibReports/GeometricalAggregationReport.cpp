#include "GeometricalAggregationReport.h"
#include <string>
#include <iostream>
#include <vector>
#include <Windows.h> // TODO this is windows only which is not what I want ultimately, have not succedded in using boost yet.
#include <StringManip.h>
#include <LogFile.h>
#include <CSVFile.h>
#include <map>

using namespace std;


void GeometricalAggregationReport::printReport(string datadir, string filename)
{
	LogFile::getInstance()->writeMessage("Starting Aggregation report for " + datadir);
	vector<string> files = getFilesWithinFolder(datadir);
	map<string, vector<double>> probabilities;
	for (unsigned int i = 0; i < files.size(); ++i)
	{
		string file = datadir + files[i];
		LogFile::getInstance()->writeMessage(files[i]);
		CSVFile csv(file);
		for (unsigned int i = 1; i < csv.fileVector.size(); ++i) // TODO a bit hard coded, skipping the header and we know we want 1,3,4,5,6 then 7 is the value
		{
			string key = csv.fileVector[i][1] + ":"; //ALA
			key += csv.fileVector[i][3] + ":"; //SS
			key += csv.fileVector[i][4] + ":"; // ANGLE
			key += csv.fileVector[i][5] + ":"; // METHOD
			key += csv.fileVector[i][6]; // chain			
			double value = atof((csv.fileVector[i][7]).c_str());
			map<string, vector<double>>::iterator iter = probabilities.find(key);
			if (iter == probabilities.end())
			{
				vector<double> oneval;
				oneval.push_back(value);
				probabilities.insert(pair < string, vector<double>>(key, oneval));
			}
			else
			{
				probabilities[key].push_back(value);
			}			
		}
	}

	//now print out the probability distribution file
	stringstream report;
	report << "ProabilityID,Values\n";
	
	for (map<string, vector<double>>::iterator iter = probabilities.begin(); iter != probabilities.end(); ++iter)
	{
		string key = iter->first;
		vector<double> vals = iter->second;
		report << key;
		for (unsigned int i = 0; i < vals.size(); ++i)
		{
			report <<"," << vals[i];
		}
		report << "\n";
	}

	ofstream outfile(filename);
	if (outfile.is_open())
	{
		outfile << report.str();
	}

}

vector<string> GeometricalAggregationReport::getFilesWithinFolder(string folder)
{//https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
	vector<string> names;
	string search_path = folder + "*.*";
	wstring wrpath = StringManip::utf8ToUtf16(search_path);
	LPCWSTR lpr = wrpath.c_str();
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(lpr, &fd);
	if (hFind != INVALID_HANDLE_VALUE) 
	{
		do {
			// read all (real) files in current folder
			// , delete '!' read other 2 default folder . and ..
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) 
			{
				wstring fs = fd.cFileName;				
				string s = StringManip::ws2s(fs);
				names.push_back(s);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}


