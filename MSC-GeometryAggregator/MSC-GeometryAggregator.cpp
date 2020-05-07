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
	string geopath = "F:\\PSUA\\ProjectData\\GeoTest4\\GeoValues\\";
	string outpath = "F:\\PSUA\\ProjectData\\GeoTest4\\Aggregation\\";
	string logfile = "F:\\PSUA\\ProjectData\\GeoTest4\\agglogger.txt";
	string pdblist = "F:\\PSUA\\Code\\PSU-ALPHA\\Project\\PDBData\\CandidateList1.csv";
	
	bool success = LogFile::getInstance()->setLogFile(logfile, outpath);
	LogFile::getInstance()->writeMessage("********** Starting Geometry Aggregation for PSU-BETA **************");

	CSVFile pdblistfile(pdblist, ",", true);	

	vector<string> datafiles;
	for (unsigned int i = 1; i < pdblistfile.fileVector.size(); ++i)
	{
		string pdb = pdblistfile.fileVector[i][0];
		string filename = geopath + pdb + "_A_geo.txt";
		//string filename2 = geopath + "Dih" + pdb + "_A_geo.txt"; //tmp code as missed off phi psi and omega
		datafiles.push_back(filename);
		//datafiles.push_back(filename2);
	}

	/*
	in this first version I am only going to load up the A version of the occupants would I deal with them all sperately or all together
	*/	
	LogFile::getInstance()->writeMessage("Outputting geometric features to, outfile=" + geopath);
	GeometricalAggregationReport gar;
	/*Memory versus time
	My laptop can't cope with all the data in internal memory at once
	So I am going through every file and pulling out the data for each amino acid in turn
	It takes longer but releases memory
	*/
	vector<string> aminos;	
	aminos.push_back("ALA");
	aminos.push_back("CYS");
	aminos.push_back("ASP");
	aminos.push_back("GLU");
	aminos.push_back("PHE");
	aminos.push_back("GLY");
	aminos.push_back("HIS");
	aminos.push_back("ILE");
	aminos.push_back("LYS");
	aminos.push_back("LEU");
	aminos.push_back("MET");
	aminos.push_back("ASN");
	aminos.push_back("PRO");
	aminos.push_back("GLN");
	aminos.push_back("ARG");
	aminos.push_back("SER");
	aminos.push_back("THR");
	aminos.push_back("VAL");
	aminos.push_back("TRP");
	aminos.push_back("TYR");
	for (unsigned int a = 0; a < aminos.size(); ++a)
	{
		LogFile::getInstance()->writeMessage("-- Aggregating " + aminos[a] + "--");
		gar.printReport(datafiles, outpath, aminos[a]);
	}


	cout << "finished";
}

