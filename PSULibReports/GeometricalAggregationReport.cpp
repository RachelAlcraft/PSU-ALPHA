#include "GeometricalAggregationReport.h"
#include <string>
#include <iostream>
#include <vector>
#include <Windows.h> // TODO this is windows only which is not what I want ultimately, have not succedded in using boost yet.

#include <LogFile.h>
#include <CSVFile.h>
#include <map>
#include <FoldersFiles.h>

using namespace std;


void GeometricalAggregationReport::printReport(string datadir)
{
	LogFile::getInstance()->writeMessage("Starting Aggregation report for " + datadir);
	vector<string> files = FoldersFiles::getFilesWithinFolder(datadir);
	map<string,map<string, vector<double>>> probabilities;
	for (unsigned int i = 0; i < files.size(); ++i)
	{
		stringstream status;
		status << i << " out of " << files.size() << " ";
		string file = datadir + files[i];
		LogFile::getInstance()->writeMessage(status.str() + files[i]);
		try
		{					
			CSVFile csv(file);
			for (unsigned int i = 1; i < csv.fileVector.size(); ++i) // TODO a bit hard coded, skipping the header and we know we want 1,3,4,5,6 then 7 is the value
			{
				if (csv.fileVector[i].size() > 7)
				{
					string key = csv.fileVector[i][1] + ":"; //ALA
					key += csv.fileVector[i][3] + ":"; //SS
					key += csv.fileVector[i][4] + ":"; // ANGLE
					key += csv.fileVector[i][5] + ":"; // METHOD
					key += csv.fileVector[i][6]; // chain			
					double value = atof((csv.fileVector[i][7]).c_str());
					map<string, map<string, vector<double>>>::iterator aiter = probabilities.find(csv.fileVector[i][1]);
					if (aiter == probabilities.end())
					{
						map<string, vector<double>> emptydata;
						probabilities.insert(pair<string, map<string, vector<double>>>(csv.fileVector[i][1], emptydata));
					}
					//Now it is there for sure
					map<string, map<string, vector<double>>>::iterator biter = probabilities.find(csv.fileVector[i][1]);
					map<string, vector<double>>::iterator citer = biter->second.find(key);
					if (citer == biter->second.end())
					{
						vector<double> emptydata;
						probabilities[csv.fileVector[i][1]].insert(pair<string, vector<double>>(key, emptydata));
					}
					//Now it is there for sure
					probabilities[csv.fileVector[i][1]][key].push_back(value);
				}
				else
				{
					LogFile::getInstance()->writeMessage("Invalid data for creating geometric aggregation");
				}
			}
		}
		catch (...)
		{
			LogFile::getInstance()->writeMessage("!!!Error");
		}
	}

	//now print out the probability distribution file per amino acid so it is easy to look at
	
	
	for (map<string,map<string, vector<double>>>::iterator iter = probabilities.begin(); iter != probabilities.end(); ++iter)
	{
		try
		{
			//TODO better to make it a dataframe by having
			// AA,Exmp_meth,geo type, SS,atoms,values as a : delim list
			// Then the prob dist can be got at whatever preferre granularity
			stringstream report;
			report << "ProabilityID,Values\n";

			string aa = iter->first;
			for (map<string, vector<double>>::iterator biter = iter->second.begin(); biter != iter->second.end(); ++biter)
			{
				string key = biter->first;
				vector<double> vals = biter->second;
				report << key;
				for (unsigned int i = 0; i < vals.size(); ++i)
				{
					report << "," << vals[i];
				}
				report << "\n";
			}

			string filename = datadir + aa + "_geoprobdist.csv";
			ofstream outfile(filename);
			if (outfile.is_open())
			{
				LogFile::getInstance()->writeMessage("Printing " + filename);
				outfile << report.str();
			}
		}
		catch (...)
		{
			LogFile::getInstance()->writeMessage("!!! Error printing aggregation");
		}
	}

	

}

/*vector<string> GeometricalAggregationReport::getFilesWithinFolder(string folder)
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
}*/


