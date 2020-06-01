// MSC-CreatePDBGeometry.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <PDBFile.h>
#include <ProteinManager.h>
#include <GeometricalDataReport.h>
#include <CSVFile.h>
#include <LogFile.h>
#include <CAlphaReport.h>

using namespace std;

int main()
{

	// ***********************************************
	// * USER UNPUT ******** 
	string pdblist = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_highres_90_annotated.csv";
	//string pdblist = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2019_95_annotated.csv";
	
	string pdbdir = "F:\\Code\\BbkTransfer\\pdbfiles\\pdbdata\\";
	string outputdir = "F:\\Code\\BbkDatabase\\tbl2_geo_contact\\DataSets\\Version1\\29May20\\";

	string contactA = "N";
	string contactB = "O";

	double max_distance = 8;




	//***************************************************************************
   // Execution begins

	//string geoOutput = outputdir + "GeoValues\\";
	//string atomOutput = outputdir + "Atoms\\";
	//string calclist = outputdir + "GeoCalc.csv";
	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";

	bool success = LogFile::getInstance()->setLogFile("logger.txt", "");
	LogFile::getInstance()->writeMessage("********** Starting Geometry calculations for PSU-BETA **************");

	ProteinManager::getInstance()->createAminoAcidData(aadata);


	CSVFile pdblistfile(pdblist, ",", true);

	CAlphaReport car;
	
	for (unsigned int i = 1; i < pdblistfile.fileVector.size(); ++i)
	{		
		string pdb = pdblistfile.fileVector[i][0];
		//if (pdb == "6SUP")
		{
			stringstream ss;
			ss << "-----" << i << " -" << pdb << "--------";
			LogFile::getInstance()->writeMessage(ss.str());

			//string pdb = "3BVX";
			PDBFile* pdbf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, pdbdir + pdb + ".pdb");
			pdbf->loadData();
			pdbf->loadAtoms();
			pdbf->loadBonds();

			car.printReport(pdbf, "A", "", "", contactA, contactB, outputdir + pdb + "_" + contactA + "_" + contactB + "_contact.csv",max_distance);

			ProteinManager::getInstance()->deletePdbs();//keep memory clear
		}
	}
	cout << "Finished";
}


