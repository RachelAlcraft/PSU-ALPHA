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
	
	string pdbdir = "F:\\PSUA\\ProjectData\\HighResFiles\\";
	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";
	
	string logpath = "F:\\PSUA\\ProjectData\\GeoTest5\\";
	string logfile = "F:\\PSUA\\ProjectData\\GeoTest5\\logger.txt";
	string geoOutput = "F:\\PSUA\\ProjectData\\GeoTest5\\GeoValues\\";
	string pdblist = "F:\\PSUA\\ProjectData\\GeoTest5\\CandidateList1.csv";
	string geo = "F:\\PSUA\\ProjectData\\GeoTest5\\GeoFeatures.csv";
	
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


