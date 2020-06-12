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
	vector<string> pdblists;
	pdblists.push_back("F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_highres_90_annotated.csv");
	pdblists.push_back("F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\list_2019_95_annotated.csv");
	
	
	string pdbdir = "F:\\Code\\BbkTransfer\\pdbfiles\\pdbdata\\";
	string outputdir = "F:\\Code\\BbkDatabase\\tbl2_geo_contact\\DataSets\\Version1\\04Jun20\\";

	vector<string> contactsA;
	vector<string> contactsB;
	
	contactsA.push_back("CA");
	contactsB.push_back("CA");

	contactsA.push_back("CB");
	contactsB.push_back("CB");

	double max_distance = 6.1;




	//***************************************************************************
   // Execution begins

	//string geoOutput = outputdir + "GeoValues\\";
	//string atomOutput = outputdir + "Atoms\\";
	//string calclist = outputdir + "GeoCalc.csv";
	string aadata = "F:\\PSUA\\Code\\PSU-ALPHA\\Config\\data_aminoinfo.csv";

	bool success = LogFile::getInstance()->setLogFile("logger.txt", "");
	LogFile::getInstance()->writeMessage("********** Starting Geometry calculations for PSU-BETA **************");

	ProteinManager::getInstance()->createAminoAcidData(aadata);
	
	for (unsigned int p = 0; p < pdblists.size(); ++p)
	{
		string pdblist = pdblists[p];

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
				try
				{
					PDBFile* pdbf = ProteinManager::getInstance()->getOrAddPDBFile(pdb, pdbdir + pdb + ".pdb");
					pdbf->loadData();
					pdbf->loadAtoms();
					pdbf->loadBonds();

					for (unsigned int c = 0; c < contactsA.size(); ++c)
					{
						string contactA = contactsA[c];
						string contactB = contactsB[c];

						LogFile::getInstance()->writeMessage("---------" + contactA + "-" + contactB);

						car.printReport(pdbf, "A", "", "", contactA, contactB, outputdir + pdb + "_" + contactA + "_" + contactB + "_contact.csv", max_distance);
					}
					ProteinManager::getInstance()->deletePdbs();//keep memory clear
				}
				catch (...)
				{
					LogFile::getInstance()->writeMessage("!!! There was an error with this file");
				}
			}
		}
	}
	cout << "Finished";
}


