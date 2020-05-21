// MSC-CreatePDBGeometry.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <PDBFile.h>
#include <ProteinManager.h>
#include <GeometricalDataReport.h>
#include <CSVFile.h>
#include <LogFile.h>

using namespace std;

int main()
{
	string setname = "\\ProjectData\\GeoPhi\\";

	string logpath = "F:\\PSUA" + setname;
	string logfile = "F:\\PSUA" + setname + "logger.txt";
	string geoOutput = "F:\\PSUA" + setname + "GeoValues\\";
	string pdblist = "F:\\PSUA" + setname + "CandidateList1.csv";
	string geo = "F:\\PSUA" + setname + "GeoCalcPhi.csv";

	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";
	string pdbdir = "F:\\PSUA\\ProjectData\\HighResFiles\\";
	
	bool success = LogFile::getInstance()->setLogFile(logfile, logpath);
	LogFile::getInstance()->writeMessage("********** Starting Geometry calculations for PSU-BETA **************");
	
	ProteinManager::getInstance()->createAminoAcidData(aadata);

	
	CSVFile pdblistfile(pdblist, ",", true);	
	CSVFile geoFile(geo, ",", true);
	
	GeometricalDataReport gdr(geoFile.fileVector);

	for (unsigned int i = 1; i < pdblistfile.fileVector.size(); ++i)
	{
		string pdb = pdblistfile.fileVector[i][0];		
		stringstream ss;
		ss << "-----" << i << " -" << pdb << "--------";
		LogFile::getInstance()->writeMessage(ss.str());

		//string pdb = "3BVX";
		PDBFile* pdbf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, pdbdir + pdb + ".pdb");
		pdbf->loadData();
		pdbf->loadAtoms();
		pdbf->loadBonds();		
		gdr.printReport(pdbf, geoOutput + pdb, "", "");

		ProteinManager::getInstance()->deletePdbs();//keep memory clear
	}
	cout << "Finished";    
}


