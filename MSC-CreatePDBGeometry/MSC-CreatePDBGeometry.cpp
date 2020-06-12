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

	// ***********************************************
	// * USER UNPUT ******** 
	vector<string> pdblists;
	//pdblists.push_back("F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_highres_90_annotated.csv");
	pdblists.push_back("F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2019_95_annotated.csv");
	//pdblists.push_back("list_test1.csv");
	//pdblists.push_back("list_test2.csv");
	
	string pdbdir = "F:\\Code\\BbkTransfer\\pdbfiles\\pdbdata\\";		
	string outputdir = "F:\\Code\\BbkDatabase\\tbl2_geo_measure_and_atom\\DataSets\\Version4\\04June20_2\\";

	bool runCore = true;
	bool runExtra = false;
	bool runAtoms = false;

	
	

	 //***************************************************************************
	// Execution begins

	//string geoOutput = outputdir + "GeoValues\\";
	//string atomOutput = outputdir + "Atoms\\";
	string calclist = outputdir + "GeoCalc.csv";		
	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";
	
	bool success = LogFile::getInstance()->setLogFile("logger.txt", "");
	LogFile::getInstance()->writeMessage("********** Starting Geometry calculations for PSU-BETA **************");
	
	ProteinManager::getInstance()->createAminoAcidData(aadata);
	
	CSVFile geoFileCore(calclist, ",", true);
	
	for (unsigned int p = 0; p < pdblists.size(); ++p)
	{
		string pdblist = pdblists[p];
		CSVFile pdblistfile(pdblist, ",", true);


		GeometricalDataReport gdr(geoFileCore.fileVector);
		

		for (unsigned int i = 1; i < pdblistfile.fileVector.size(); ++i)
		{
			string pdb = pdblistfile.fileVector[i][0];
			//if(pdb == "6V98")
			{
				stringstream ss;
				ss << "-----" << i << " -" << pdb << "--------";
				LogFile::getInstance()->writeMessage(ss.str());

				//string pdb = "3BVX";
				PDBFile* pdbf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, pdbdir + pdb + ".pdb");
				pdbf->loadData();
				pdbf->loadAtoms();
				pdbf->loadBonds();

				gdr.printReport(pdbf, outputdir, runCore, runExtra, runAtoms);
				ProteinManager::getInstance()->deletePdbs();//keep memory clear
			}
		}
	}
	cout << "Finished";    
}


