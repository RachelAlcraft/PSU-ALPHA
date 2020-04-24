// MSC-GeometryAggregator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <PDBFile.h>
#include <ProteinManager.h>
#include <GeometricalDataReport.h>
#include <CSVFile.h>
#include <LogFile.h>
#include <GeometricalAggregationReport.h>

int main()
{
	string geopath = "F:\\PSUA\\ProjectData\\GeoData\\";
	string outpath = "F:\\PSUA\\ProjectData\\Aggregation\\";
	string logfile = "F:\\PSUA\\ProjectData\\Aggregation\\logger.txt";
	string pdblist = "F:\\PSUA\\Code\\PSU-ALPHA\\Project\\PDBData\\CandidateList1.csv";
	
	bool success = LogFile::getInstance()->setLogFile(logfile, outpath);
	LogFile::getInstance()->writeMessage("********** Starting Geometry Aggregation for PSU-BETA **************");

	CSVFile pdblistfile(pdblist, ",", true);	

	vector<string> datafiles;
	for (unsigned int i = 1; i < pdblistfile.fileVector.size(); ++i)
	{
		string pdb = pdblistfile.fileVector[i][0];
		string filename = geopath + pdb + "_A_geo.txt";
		datafiles.push_back(filename);
	}

	/*
	in this first version I am only going to load up the A version of the occupants would I deal with them all sperately or all together
	*/	
	LogFile::getInstance()->writeMessage("Outputting geometric features to, outfile=" + geopath);
	GeometricalAggregationReport gar;
	gar.printReport(datafiles,outpath);	


	cout << "finished";
}

